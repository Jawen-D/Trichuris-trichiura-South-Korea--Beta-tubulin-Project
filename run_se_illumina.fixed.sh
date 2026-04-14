#!/usr/bin/env bash
set -euo pipefail

# ---------------------------
# Single-end Illumina pipeline (bcftools + SnpEff)
# - Aligns with BWA-MEM (adds read group so the sample name is correct)
# - Calls variants with bcftools (adds FORMAT/DP and AD)
# - Filters, normalizes, renames contig to AF034219.1 (for SnpEff match)
# - Annotates with SnpEff (logs separated, valid VCF)
# - Exports a tidy TSV
# ---------------------------

usage() {
  echo "Usage: $0 -s SAMPLE -q reads.fastq[.gz] [-p THREADS (default 8)]"
  exit 1
}

# ---- config (run from your project root where ref.fasta lives) ----
PROJ="$(pwd)"
REF="${PROJ}/ref.fasta"
SNPEFF_DIR="/opt/anaconda3/envs/varcall/share/snpeff-5.2-1"
GENOME="Trichuris_custom"
THREADS=8

# ---- args ----
S="ER54"; Q="sample/ER54/ER54_rv.fastq"
while getopts ":s:q:p:" opt; do
  case $opt in
    s) S="$OPTARG" ;;
    q) Q="$OPTARG" ;;
    p) THREADS="$OPTARG" ;;
    *) usage ;;
  esac
done
[[ -z "${S}" || -z "${Q}" ]] && usage

# ---- sanity checks ----
[[ -f "${REF}" ]] || { echo "ERROR: Missing ${REF}"; exit 2; }
[[ -f "${Q}"  ]] || { echo "ERROR: FASTQ not found: ${Q}"; exit 2; }
command -v bwa >/dev/null || { echo "ERROR: bwa not in PATH"; exit 2; }
command -v samtools >/dev/null || { echo "ERROR: samtools not in PATH"; exit 2; }
command -v bcftools >/dev/null || { echo "ERROR: bcftools not in PATH"; exit 2; }
[[ -f "${SNPEFF_DIR}/snpEff.jar" ]] || { echo "ERROR: snpEff.jar not found in ${SNPEFF_DIR}"; exit 2; }

mkdir -p "${PROJ}/bam" "${PROJ}/vcf"

echo ">>> SAMPLE=${S}"
echo ">>> FASTQ =${Q}"
echo ">>> REF   =${REF}"
echo

# ---- 1) reference indexes (idempotent) ----
samtools faidx "${REF}"
# Make BWA index if missing
[[ -f "${REF}.bwt" ]] || bwa index "${REF}"

# ---- 2) align (add read group so the VCF sample is '${S}') ----
BAM="${PROJ}/bam/${S}.sorted.bam"
echo ">> Aligning with BWA-MEM..."
bwa mem -t "${THREADS}" -R "@RG\tID:${S}\tSM:${S}" "${REF}" "${Q}" \
  | samtools sort -@ "${THREADS}" -o "${BAM}"
samtools index "${BAM}"

# ---- 3) variant calling (ensure AD/DP present) ----
RAW="${PROJ}/vcf/${S}.raw.vcf.gz"
echo ">> Calling variants with bcftools..."
bcftools mpileup -f "${REF}" -q 20 -Q 20 -a FORMAT/DP,AD "${BAM}" \
| bcftools call -mv -Oz -o "${RAW}"
bcftools index -f "${RAW}"

# ---- 4) filter + normalize ----
FILT="${PROJ}/vcf/${S}.filtered.vcf.gz"
FN="${PROJ}/vcf/${S}.filtered.norm.vcf.gz"
echo ">> Filtering (QUAL>=30 && DP>=10) and normalizing..."
bcftools +fill-tags "${RAW}" -- -t AF,AC,AN \
| bcftools filter -i 'QUAL>=30 && FORMAT/DP>=10' -Oz -o "${FILT}"
bcftools index -f "${FILT}"
bcftools norm -f "${REF}" -m -both "${FILT}" -Oz -o "${FN}"
bcftools index -f "${FN}"

# ---- 5) contig rename to match SnpEff DB (AF034219 -> AF034219.1) ----
REN="${PROJ}/vcf/${S}.filtered.norm.renamed.vcf.gz"
MAP="${PROJ}/vcf/chr_rename.txt"
printf "AF034219\tAF034219.1\nAF034219.1\tAF034219.1\n" > "${MAP}"
echo ">> Ensuring contig name matches SnpEff DB (AF034219.1)..."
# Try rename; if it fails (already matching), just copy
if ! bcftools annotate --rename-chrs "${MAP}" -Oz -o "${REN}" "${FN}" ; then
  cp "${FN}" "${REN}"
fi
bcftools index -f "${REN}"

# ---- 6) SnpEff annotation (keep logs out of the VCF) ----
ANN_PLAIN="${PROJ}/vcf/${S}.annotated.vcf"
ANN="${PROJ}/vcf/${S}.annotated.vcf.gz"
LOG="${PROJ}/vcf/${S}.snpeff.log"
echo ">> Annotating with SnpEff (${GENOME})..."
set -o pipefail
java -Xmx4g -jar "${SNPEFF_DIR}/snpEff.jar" -v "${GENOME}" "${REN}" \
  1> "${ANN_PLAIN}" 2> "${LOG}"
bgzip -f "${ANN_PLAIN}"
bcftools index -f "${ANN}"

# Sanity: ANN header must exist
if ! bcftools view -h "${ANN}" | grep -q '^##INFO=<ID=ANN'; then
  echo "ERROR: Annotated VCF lacks ANN header. See ${LOG}" >&2
  exit 3
fi

# ---- 7) tidy TSV (ANN parsed) ----
TSV="${PROJ}/vcf/${S}.annotated.parsed.tsv"
echo ">> Writing parsed table: ${TSV}"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/ANN\n' "${ANN}" \
| awk -F'\t' 'BEGIN{OFS="\t"}{
  # Take first ANN record if multiple
  split($6,all,","); split(all[1],ann,"|");
  effect=(length(ann)>=2?ann[2]:"");
  impact=(length(ann)>=3?ann[3]:"");
  gene=(length(ann)>=4?ann[4]:"");
  tx=(length(ann)>=7?ann[7]:"");
  aa=(length(ann)>=11?ann[11]:"");
  print $1,$2,$3,$4,$5,effect,impact,gene,tx,aa
}' > "${TSV}"

echo
echo "Done."
echo "BAM:         ${BAM}"
echo "Raw VCF:     ${RAW}"
echo "Filtered:    ${FILT}"
echo "Norm:        ${FN}"
echo "Renamed:     ${REN}"
echo "Annotated:   ${ANN}"
echo "Parsed TSV:  ${TSV}"
echo "SnpEff log:  ${LOG}"
