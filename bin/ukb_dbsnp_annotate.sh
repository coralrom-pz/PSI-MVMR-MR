
```bash
#!/usr/bin/env bash
# =========================================================================================
# SCRIPT — dbSNP RSID Appender for UKB-PPP per-chromosome pQTL tables
# -----------------------------------------------------------------------------------------
# Purpose:
#   Append RSIDs (from dbSNP VCF) to UKB per-chrom discovery tables inside each protein dir.
#
# Usage:
#   chmod +x scripts/ukb_ppp_dbsnp_rsid.sh
#   UKB_ROOT="project/data/UKB-PPP" \
#   DBSNP="project/resources/dbsnp/dbsnp.vcf.gz" \
#   THREADS=4 \
#   scripts/ukb_ppp_dbsnp_rsid.sh
#
# Inputs per protein folder (e.g., IL1B_P01584_OID20427_v1_Inflammation/):
#   discovery_chr1_*.txt(.gz), discovery_chr2_*.txt(.gz), ...
#
# Outputs per protein folder:
#   processed/IL1B_chr1_Inflammation_with_rsid.txt
#   IL1B_with_rsid_allchr.txt   (merged with single header)
#
# Requirements:
#   bcftools, awk, sed, mktemp, xargs, gunzip
# =========================================================================================
set -euo pipefail

: "${UKB_ROOT:="project/data/UKB-PPP"}"
: "${DBSNP:="project/resources/dbsnp/dbsnp.vcf.gz"}"
: "${THREADS:=4}"

need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found."; exit 1; }; }
need bcftools; need awk; need sed; need mktemp; need xargs; need gunzip

[[ -d "$UKB_ROOT" ]] || { echo "ERROR: UKB_ROOT not found: $UKB_ROOT"; exit 1; }
[[ -f "$DBSNP" && -f "${DBSNP}.tbi" ]] || { echo "ERROR: dbSNP VCF or index missing: $DBSNP(.tbi)"; exit 1; }

echo "Using UKB_ROOT = $UKB_ROOT"
echo "Using DBSNP    = $DBSNP"

# ---- Build contig name map from dbSNP header (supports '1', 'chr1', 'NC_000001.11') ----
CHRMAP="$(mktemp -t dbsnp_chrmap_XXXXXX.tsv)"
bcftools view -h "$DBSNP" | awk -F'[<>=,]' '
  BEGIN{OFS="\t"}
  /^##contig=/ {
    id="";
    for(i=1;i<=NF;i++) if($i=="ID") { id=$(i+1) }
    if(id ~ /^NC_[0-9]+[.][0-9]+$/) {
      num=id; sub(/^NC_0+/, "", num); sub(/[.].*$/, "", num)
      if(num=="23") key="X"; else if(num=="24") key="Y"; else key=num
      print key, id
    } else if (id ~ /^chr[0-9XY]+$/) {
      key=id; sub(/^chr/, "", key); print key, id
    } else if (id ~ /^[0-9XY]+$/) {
      print id, id
    }
  }' > "$CHRMAP"

if ! grep -qE '^[0-9XY]+\t' "$CHRMAP"; then
  echo "ERROR: Could not infer dbSNP contig naming from header. Inspect: $CHRMAP"
  exit 1
fi
echo "Detected contig naming (first few):"; head -n 5 "$CHRMAP"

annotate_table_to_out() {
  local orig="$1"   # temp path to the original table (header present)
  local out="$2"    # destination .txt with RSID appended

  local tmpdir; tmpdir="$(mktemp -d -t ukb_rsid_XXXXXX)"
  trap 'rm -rf "$tmpdir"' RETURN
  local vcf="${tmpdir}/min.vcf"
  local annvcf="${tmpdir}/min.rsid.vcf"
  local mapfile_tsv="${tmpdir}/pos2rsid.tsv"

  # 1) Build mini VCF from original table row IDs (CHROM:POS:REF:ALT expected in col 3)
  awk 'BEGIN{OFS="\t"} NR>1{
         split($3,a,":"); if(length(a)<4) next;
         print a[1], a[2], ".", a[3], a[4], ".", "PASS", "."
       }' "$orig" \
  | { printf "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"; cat; } > "$vcf"

  # Sort the mini VCF
  bcftools sort -T "$tmpdir" -O v -o "${vcf}.sorted" "$vcf"
  mv "${vcf}.sorted" "$vcf"

  # ---- Position-only RSID lookup via bcftools query ----
  local regions_ukb="${tmpdir}/regions.ukb.tsv"
  awk 'NR>1{print $1 "\t" $2 "\t" $2}' "$orig" > "$regions_ukb"

  # If no regions, emit empty RSID col and return
  if [[ ! -s "$regions_ukb" ]]; then
    awk 'NR==1{print $0 "\tRSID"; next} {print $0 "\t"}' "$orig" > "$out"
    return 0
  fi

  # Build reverse contig map: dbSNP contig -> UKB contig (1..22/X/Y)
  local revmap="${tmpdir}/rev.tsv"
  awk 'BEGIN{OFS="\t"}{print $2,$1}' "$CHRMAP" > "$revmap"

  # Rename regions to dbSNP contigs
  local regions_dbsnp="${tmpdir}/regions.dbsnp.tsv"
  awk -v OFS="\t" '
    BEGIN{
      while((getline < ARGV[1])>0){ f[$1]=$2 } close(ARGV[1]); ARGV[1]=""
    }
    { chr=$1; if(chr in f) $1=f[chr]; print $0 }
  ' "$CHRMAP" "$regions_ukb" > "$regions_dbsnp"

  # Query dbSNP for rsIDs at those positions
  local pos2id_dbsnp="${tmpdir}/pos2id_dbsnp.tsv"
  if ! bcftools query -R "$regions_dbsnp" -f '%CHROM\t%POS\t%ID\n' "$DBSNP" > "$pos2id_dbsnp" 2> "${tmpdir}/bcf.err"; then
    echo "bcftools query failed; see ${tmpdir}/bcf.err"
    awk 'NR==1{print $0 "\tRSID"; next} {print $0 "\t"}' "$orig" > "$out"
    return 0
  fi

  # Convert dbSNP contigs back to UKB style (so we can join on UKB CHROM:POS)
  local pos2id_ukb="${tmpdir}/pos2id_ukb.tsv"
  awk -v OFS="\t" '
    BEGIN{
      while((getline < ARGV[1])>0){ r[$1]=$2 } close(ARGV[1]); ARGV[1]=""
    }
    { chr=$1; if(chr in r) chr=r[chr]; print chr, $2, $3 }
  ' "$revmap" "$pos2id_dbsnp" > "$pos2id_ukb"

  # Build a map CHROM POS -> rsID (take the first rsID if multiple exist)
  local mapfile_tsv="${tmpdir}/pos2rsid.tsv"
  awk -v OFS="\t" '{k=$1 FS $2; ids[k]=(k in ids?ids[k]";"$3:$3)} END{for(k in ids){split(k,a,"\t"); print a[1],a[2],ids[k]}}' "$pos2id_ukb" > "$mapfile_tsv"
  
  # Append RSID to original table by joining on CHROM:POS from the ID field
  awk -v OFS="\t" '
    FNR==NR { m[$1 FS $2]=$3; next }   # map: CHROM POS -> rsID
    FNR==1  { print $0, "RSID"; next } # header of the pQTL table
             {
               key = $1 FS $2;         # CHROM and GENPOS from columns 1 and 2
               val = (key in m ? m[key] : "")
               print $0, val
             }' "$mapfile_tsv" "$orig" > "$out"
}

export -f annotate_table_to_out
export DBSNP CHRMAP

process_one_file() {
  local infile="$1"   # absolute path to discovery_chr*.txt(.gz)
  local pdir="$2"     # protein directory (absolute)
  local prot_short="$3"
  local panel="$4"

  local processed_dir="${pdir}/processed"
  mkdir -p "$processed_dir"

  local base; base="$(basename "$infile")"
  local chrom; chrom="$(echo "$base" | sed -nE 's/^discovery_chr([0-9XY]+)_.*/\1/p')"
  [[ -n "${chrom:-}" ]] || chrom="NA"

  local out="${processed_dir}/${prot_short}_chr${chrom}_${panel}_with_rsid.txt"

  # copy input to temp
  local tmp; tmp="$(mktemp -t ukb_orig_XXXXXX.txt)"
  case "$infile" in
    *.gz) gunzip -c "$infile" > "$tmp" ;;
    *)    cp "$infile" "$tmp" ;;
  esac

  annotate_table_to_out "$tmp" "$out"
  rm -f "$tmp"
  echo "✔ ${prot_short} chr${chrom}: $out"
}

export -f process_one_file

# Walk each protein directory
for pdir in "$UKB_ROOT"/*/ ; do
  [[ -d "$pdir" ]] || continue
  prot_dirname="$(basename "$pdir")"                          # e.g., IL1B_P01584_OID20427_v1_Inflammation
  prot_short="$(echo "$prot_dirname" | awk -F_ '{print $1}')" # IL1B
  panel="${prot_dirname##*_}"                                 # Inflammation / Neurology / Cardiometabolic

  # Find per-chrom files; process in parallel
  file_glob_count=$(find "$pdir" -maxdepth 1 -type f \( -name "discovery_chr*.txt" -o -name "discovery_chr*.txt.gz" \) | wc -l | tr -d ' ')
  if [[ "$file_glob_count" -eq 0 ]]; then
    echo "No discovery_chr*.txt files in $pdir — skipping."
    continue
  fi

  find "$pdir" -maxdepth 1 -type f \( -name "discovery_chr*.txt" -o -name "discovery_chr*.txt.gz" \) -print0 \
  | xargs -0 -I{} -P "$THREADS" bash -c 'process_one_file "$1" "$2" "$3" "$4"' _ "{}" "$pdir" "$prot_short" "$panel"

  # Merge per-chrom outputs for this protein (single header)
  processed_dir="${pdir}/processed"
  merged="${pdir}/${prot_short}_with_rsid_allchr.txt"
  if ls "${processed_dir}/${prot_short}_chr"*"_${panel}_with_rsid.txt" >/dev/null 2>&1; then
    ls -1 "${processed_dir}/${prot_short}_chr"*"_${panel}_with_rsid.txt" \
      | sed -E 's/_chrX_/_chr23_/; s/_chrY_/_chr24_/' \
      | sort \
      | sed -E 's/_chr23_/_chrX_/; s/_chr24_/_chrY_/' \
      | {
          read first
          head -n 1 "$first" > "$merged"
          echo "$first"
          cat
        } | while read -r f; do tail -n +2 "$f" >> "$merged"; done
    echo "Merged → $merged"
  else
    echo "No processed files to merge for $prot_dirname"
  fi
done

echo "All done."
