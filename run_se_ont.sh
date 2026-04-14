#!/usr/bin/env bash
set -euo pipefail

usage(){ echo "Usage: $0 -s SAMPLE -q reads.fastq[.gz] [-p THREADS]"; exit 1; }

# config (edit if your paths differ)
PROJ="$(pwd)"
REF="${PROJ}/ref.fasta"
SNPEFF_DIR="/opt/anaconda3/envs/varcall/share/snpeff-5.2-1"
GENOME="Trichuris_custom"
THREADS=8

S="TT18"; Q=""
while getopts ":s:q:p:" opt; do
  case "$opt" in
    s) S="$OPTARG" ;;
    q) Q="$OPTARG" ;;
    p) THREADS="$OPTARG" ;;
    *) usage ;;
  esac
done
[ -z "${S}" ] && usage
[ -z "${Q}" ] && usage

# tools & inputs
command -v minimap2 >/dev/null || { echo "ERROR: minimap2 not in PATH"; exit 2; }
command -v samtools  >/dev/null || { echo "ERROR: samtools not in PATH"; exit 2; }
command -v bcftools  >/dev/null || { echo "ERROR: bcftools not in PATH"; exit 2; }
[ -f "${SNPEFF_DIR}/snpEff.jar" ] || { echo "ERROR: snpEff.jar not found in ${SNPEFF_DIR}"; exit 2; }
[ -f "${REF}" ] || { echo "ERROR: Missing ${REF}"; exit 2; }
[ -f "${Q}" ]   || { echo "ERROR: FASTQ not found: ${Q}"; exit 2; }

mkdir -p "${PROJ}/bam" "${PROJ}/vcf"

# 1) reference index (idempotent)
samtools faidx "${REF}"

# 2) align (ONT) + sort/index
BAM="${PROJ}/bam/${S}.sorted.bam"
minimap2 -t "${THREADS}" -ax map-ont -R "@RG\tID:${S}\tSM:${S}" "${REF}" "${Q}" \
  | samtools sort -@ "${THREADS}" -o "${BAM}"
samtools index "${BAM}"

# 3) call variants (bcftools; keep AD/DP in FORMAT)
RAW="${PROJ}/vcf/${S}.raw.vcf.gz"
bcftools mpileup -f "${REF}" -q 10 -Q 7 -a FORMAT/DP,AD "${BAM}" \
| bcftools call -mv -Oz -o "${RAW}"
bcftools index -f "${RAW}"

# 4) filter + normalize
FILT="${PROJ}/vcf/${S}.filtered.vcf.gz"
FN="${PROJ}/vcf/${S}.filtered.norm.vcf.gz"
bcftools +fill-tags "${RAW}" -- -t AF,AC,AN \
| bcftools filter -i 'QUAL>=20 && FORMAT/DP>=15' -Oz -o "${FILT}"
bcftools index -f "${FILT}"
bcftools norm -f "${REF}" -m -both "${FILT}" -Oz -o "${FN}"
bcftools index -f "${FN}"

# 5) contig rename to match SnpEff DB (AF034219 -> AF034219.1), no-op if already matching
REN="${PROJ}/vcf/${S}.filtered.norm.renamed.vcf.gz"
MAP="${PROJ}/vcf/chr_rename.txt"
printf "AF034219\tAF034219.1\nAF034219.1\tAF034219.1\n" > "${MAP}"
if ! bcftools annotate --rename-chrs "${MAP}" -Oz -o "${REN}" "${FN}" ; then
  cp "${FN}" "${REN}"
fi
bcftools index -f "${REN}"

# 6) SnpEff annotate (clean stdout/stderr), then bgzip+index, verify ANN header
ANN_PLAIN="${PROJ}/vcf/${S}.annotated.vcf"
ANN="${PROJ}/vcf/${S}.annotated.vcf.gz"
LOG="${PROJ}/vcf/${S}.snpeff.log"
set -o pipefail
java -Xmx4g -jar "${SNPEFF_DIR}/snpEff.jar" -v "${GENOME}" "${REN}" \
  1> "${ANN_PLAIN}"  2> "${LOG}"
bgzip -f "${ANN_PLAIN}"
bcftools index -f "${ANN}"

if ! bcftools view -h "${ANN}" | grep -q '^##INFO=<ID=ANN'; then
  echo "ERROR: Annotated VCF lacks ANN header. See ${LOG}" >&2
  exit 3
fi

# 7) parsed TSV (first ANN record)
TSV="${PROJ}/vcf/${S}.annotated.parsed.tsv"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/ANN\n' "${ANN}" \
| awk -F'\t' 'BEGIN{OFS="\t"}{
  split($6,all,","); split(all[1],ann,"|");
  effect=(length(ann)>=2?ann[2]:"");
  impact=(length(ann)>=3?ann[3]:"");
  gene=(length(ann)>=4?ann[4]:"");
  tx=(length(ann)>=7?ann[7]:"");
  aa=(length(ann)>=11?ann[11]:"");
  print $1,$2,$3,$4,$5,effect,impact,gene,tx,aa
}' > "${TSV}"

echo "DONE"
echo "${BAM}"
echo "${ANN}"
echo "${TSV}"

