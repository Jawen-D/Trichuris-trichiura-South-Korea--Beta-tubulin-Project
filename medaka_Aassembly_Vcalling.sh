#!/bin/bash

# ========== USER SETTINGS ==========
READS="ER10.fastq"                 # ONT basecalled FASTQ file (amplicon reads)
REF="reference.fasta"              # Full-length β-tubulin reference
SAMPLE="ER10"                   # Output file prefix
THREADS=4                          # Number of CPU threads
MODEL="r941_min_sup_g507"          # Medaka model (adjust if needed)
# ===================================

# Step 1. Index the reference for minimap2
echo "[1] Indexing reference..."
minimap2 -d ${REF%.fasta}.mmi $REF

# Step 2. Align ONT reads to reference
echo "[2] Aligning ONT reads to reference..."
minimap2 -ax map-ont ${REF%.fasta}.mmi $READS > ${SAMPLE}.sam

# Step 3. Convert SAM to sorted BAM
echo "[3] Converting and sorting BAM..."
samtools view -bS ${SAMPLE}.sam | samtools sort -o ${SAMPLE}.sorted.bam
samtools index ${SAMPLE}.sorted.bam

# Step 4. Run Medaka Consensus
echo "[4] Running Medaka consensus to generate consensus FASTA..."
medaka_consensus \
  -i $READS \
  -d $REF \
  -o ${SAMPLE}_medaka \
  -t $THREADS \
  -m $MODEL

# Step 5. Generate VCF from Medaka consensus output
echo "[5] Running Medaka variant calling..."
medaka_variant \
  -i ${SAMPLE}_medaka/calls_to_draft.bam \
  -f $REF \
  -o ${SAMPLE}_medaka/variants \
  -t $THREADS

echo "✅ Done."
echo "Output files:"
echo "  → ${SAMPLE}_medaka/consensus.fasta"
echo "  → ${SAMPLE}_medaka/variants/round_0.vcf"
echo "  → ${SAMPLE}.sorted.bam"

