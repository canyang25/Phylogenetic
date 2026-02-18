#!/usr/bin/env bash
set -euo pipefail

# Build influenza segment consensus sequences from raw FASTQ reads.

THREADS=8
MIN_DEPTH=10
RAW_DIR="data/raw"
META_FILE="data/metadata/srr_ids.txt"
OUT_DIR="results"
REF_DIR="data/reference"
SEGMENTS=(PB2 PB1 PA HA NP NA M NS)

usage() {
  cat <<'EOF'
Usage:
  bash scripts/preprocess_consensus.sh [options]

Options:
  --threads <int>      Number of CPU threads (default: 8)
  --min_depth <int>    Minimum depth for consensus masking (default: 10)
  --raw_dir <path>     Directory containing SRR FASTQ files (default: data/raw)
  --meta_file <path>   File with one SRR ID per line (default: data/metadata/srr_ids.txt)
  --out_dir <path>     Output root directory (default: results)
  --ref_dir <path>     Reference directory (default: data/reference)
  --help               Show this help message
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --min_depth)
      MIN_DEPTH="$2"
      shift 2
      ;;
    --raw_dir)
      RAW_DIR="$2"
      shift 2
      ;;
    --meta_file)
      META_FILE="$2"
      shift 2
      ;;
    --out_dir)
      OUT_DIR="$2"
      shift 2
      ;;
    --ref_dir)
      REF_DIR="$2"
      shift 2
      ;;
    --help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ ! -f "$META_FILE" ]]; then
  echo "ERROR: Metadata file not found: $META_FILE" >&2
  exit 1
fi

for tool in curl awk sort head cut grep wc cat minimap2 samtools bcftools; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "ERROR: Required tool not found in PATH: $tool" >&2
    exit 1
  fi
done

mkdir -p "$REF_DIR" \
  "$OUT_DIR/reference_screen" \
  "$OUT_DIR/bam" \
  "$OUT_DIR/variants" \
  "$OUT_DIR/masks" \
  "$OUT_DIR/consensus" \
  "$OUT_DIR/msa/input"

echo "Checking FASTQ inputs..."
while read -r srr; do
  [[ -z "$srr" ]] && continue
  r1="${RAW_DIR}/${srr}_1.fastq"
  r2="${RAW_DIR}/${srr}_2.fastq"
  if [[ ! -f "$r1" || ! -f "$r2" ]]; then
    echo "ERROR: Missing FASTQ pair for ${srr}" >&2
    echo "Expected: $r1 and $r2" >&2
    exit 1
  fi
done < "$META_FILE"

H1_IDS="NC_026438,NC_026435,NC_026437,NC_026433,NC_026436,NC_026434,NC_026431,NC_026432"
H3_IDS="NC_007373,NC_007372,NC_007371,NC_007366,NC_007369,NC_007368,NC_007367,NC_007370"

echo "Downloading influenza references..."
curl -fsSL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${H1_IDS}&rettype=fasta&retmode=text" > "${REF_DIR}/h1_raw.fa"
curl -fsSL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${H3_IDS}&rettype=fasta&retmode=text" > "${REF_DIR}/h3_raw.fa"

echo "Normalizing segment headers..."
for infa in "${REF_DIR}/h1_raw.fa" "${REF_DIR}/h3_raw.fa"; do
  awk '
  BEGIN{
    m["NC_026438"]="PB2"; m["NC_026435"]="PB1"; m["NC_026437"]="PA"; m["NC_026433"]="HA";
    m["NC_026436"]="NP";  m["NC_026434"]="NA";  m["NC_026431"]="M";  m["NC_026432"]="NS";
    m["NC_007373"]="PB2"; m["NC_007372"]="PB1"; m["NC_007371"]="PA"; m["NC_007366"]="HA";
    m["NC_007369"]="NP";  m["NC_007368"]="NA";  m["NC_007367"]="M";  m["NC_007370"]="NS";
  }
  /^>/{acc=$1; sub(/^>/,"",acc); sub(/\..*/,"",acc); print ">" m[acc]; next}
  {print}
  ' "$infa" > "${infa/_raw.fa/_ref.fa}"
done

echo -e "sample\tref\tmapped_pct\tmapped_reads\ttotal_reads" > "${OUT_DIR}/reference_screen/mapping_summary.tsv"

echo "Reference screening: H1 vs H3..."
while read -r srr; do
  [[ -z "$srr" ]] && continue
  for r in h1 h3; do
    minimap2 -ax sr -t "$THREADS" "${REF_DIR}/${r}_ref.fa" "${RAW_DIR}/${srr}_1.fastq" "${RAW_DIR}/${srr}_2.fastq" \
      | samtools sort -@ "$THREADS" -o "${OUT_DIR}/reference_screen/${srr}.${r}.bam"
    samtools index "${OUT_DIR}/reference_screen/${srr}.${r}.bam"
    total=$(samtools view -c "${OUT_DIR}/reference_screen/${srr}.${r}.bam")
    mapped=$(samtools view -c -F 260 "${OUT_DIR}/reference_screen/${srr}.${r}.bam")
    pct=$(awk -v m="$mapped" -v t="$total" 'BEGIN{if(t==0)print 0; else printf "%.4f",100*m/t}')
    echo -e "${srr}\t${r}\t${pct}\t${mapped}\t${total}" >> "${OUT_DIR}/reference_screen/mapping_summary.tsv"
  done
done < "$META_FILE"

BEST_REF=$(
  awk 'NR>1{sum[$2]+=$3; n[$2]++} END{for(r in sum) printf "%s\t%.6f\n",r,sum[r]/n[r]}' \
    "${OUT_DIR}/reference_screen/mapping_summary.tsv" \
    | sort -k2,2nr | head -n1 | cut -f1
)
echo "Best reference set: ${BEST_REF}" | tee "${OUT_DIR}/reference_screen/best_reference.txt"

REF_FA="${REF_DIR}/${BEST_REF}_ref.fa"
SEG_REF_DIR="${OUT_DIR}/reference_screen/${BEST_REF}_segment_refs"
mkdir -p "$SEG_REF_DIR"
echo "Preparing per-segment reference FASTAs from ${REF_FA}..."
for seg in "${SEGMENTS[@]}"; do
  awk -v target="$seg" '
    /^>/{hdr=$0; sub(/^>/,"",hdr); keep=(hdr==target); if(keep) print ">" hdr; next}
    keep{print}
  ' "$REF_FA" > "${SEG_REF_DIR}/${seg}.fa"
  if ! grep -q '^>' "${SEG_REF_DIR}/${seg}.fa"; then
    echo "ERROR: Failed to extract segment ${seg} from ${REF_FA}" >&2
    exit 1
  fi
done

echo "Calling per-segment consensus with min depth ${MIN_DEPTH}..."
while read -r srr; do
  [[ -z "$srr" ]] && continue
  : > "${OUT_DIR}/consensus/${srr}.fa"
  for seg in "${SEGMENTS[@]}"; do
    seg_ref="${SEG_REF_DIR}/${seg}.fa"
    seg_bam="${OUT_DIR}/bam/${srr}.${seg}.bam"
    seg_vcf="${OUT_DIR}/variants/${srr}.${seg}.vcf.gz"
    seg_mask="${OUT_DIR}/masks/${srr}.${seg}.lowdepth.bed"
    seg_cons="${OUT_DIR}/consensus/${srr}.${seg}.fa"

    minimap2 -ax sr -t "$THREADS" "$seg_ref" "${RAW_DIR}/${srr}_1.fastq" "${RAW_DIR}/${srr}_2.fastq" \
      | samtools sort -@ "$THREADS" -o "$seg_bam"
    samtools index "$seg_bam"
    bcftools mpileup -Ou -f "$seg_ref" "$seg_bam" \
      | bcftools call --ploidy 1 -mv -Oz -o "$seg_vcf"
    bcftools index -f "$seg_vcf"
    samtools depth -aa "$seg_bam" \
      | awk -v d="$MIN_DEPTH" '$3<d{print $1"\t"$2-1"\t"$2}' > "$seg_mask"
    bcftools consensus -f "$seg_ref" -m "$seg_mask" "$seg_vcf" > "$seg_cons"

    if ! grep -q '^>' "$seg_cons"; then
      echo "ERROR: Empty consensus for sample ${srr}, segment ${seg}" >&2
      exit 1
    fi
    cat "$seg_cons" >> "${OUT_DIR}/consensus/${srr}.fa"
  done
done < "$META_FILE"

echo "Preparing segment-level FASTA inputs for MSA..."
for seg in "${SEGMENTS[@]}"; do
  : > "${OUT_DIR}/msa/input/${seg}.fa"
  while read -r srr; do
    [[ -z "$srr" ]] && continue
    awk -v sample="$srr" '/^>/{print ">"sample; next} {print}' \
      "${OUT_DIR}/consensus/${srr}.${seg}.fa" >> "${OUT_DIR}/msa/input/${seg}.fa"
  done < "$META_FILE"
  seg_n=$(grep -c '^>' "${OUT_DIR}/msa/input/${seg}.fa" || true)
  if [[ "$seg_n" -eq 0 ]]; then
    echo "ERROR: No sequences collected for segment ${seg} in ${OUT_DIR}/msa/input/${seg}.fa" >&2
    exit 1
  fi
done

echo "Preprocessing complete."
echo "Outputs:"
echo "  ${OUT_DIR}/reference_screen/mapping_summary.tsv"
echo "  ${OUT_DIR}/reference_screen/best_reference.txt"
echo "  ${OUT_DIR}/consensus/*.fa"
echo "  ${OUT_DIR}/msa/input/*.fa"
