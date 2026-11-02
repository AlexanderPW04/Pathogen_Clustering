#!/usr/bin/env bash
# OPTIONAL offline workflow: paired-end reads -> iVar variants + A/C/G/T counts.
#
# NOT run by the website, the Docker image, or the automated tests. It needs
# bwa, samtools, and ivar installed locally (for example from Bioconda) and
# enough memory and disk for raw sequencing files. Parameters follow the
# methods of Fontenele et al. (2021); confirm them against your own protocol.
#
# Usage:
#   bash offline_workflow/raw_reads_to_tables.sh SAMPLE R1.fastq.gz R2.fastq.gz \
#        reference.fasta primers.bed OUTPUT_FOLDER
set -euo pipefail

if [ "$#" -ne 6 ]; then
  echo "Usage: $0 SAMPLE R1.fastq.gz R2.fastq.gz reference.fasta primers.bed OUTPUT_FOLDER" >&2
  exit 2
fi

SAMPLE="$1"; R1="$2"; R2="$3"; REFERENCE="$4"; PRIMERS="$5"; OUT="$6"
mkdir -p "$OUT/ivar_tsv"

# 1. Map reads to the reference with BWA-MEM and sort the alignments.
[ -f "$REFERENCE.bwt" ] || bwa index "$REFERENCE"
bwa mem "$REFERENCE" "$R1" "$R2" | samtools view -b -F 4 - | samtools sort -o "$OUT/$SAMPLE.sorted.bam"
samtools index "$OUT/$SAMPLE.sorted.bam"

# 2. Soft-clip primers and quality-trim with iVar (sliding window 4, PHRED 20;
#    reads shorter than 30 nt are removed), then sort again.
ivar trim -i "$OUT/$SAMPLE.sorted.bam" -b "$PRIMERS" -p "$OUT/$SAMPLE.trimmed" -q 20 -s 4 -m 30
samtools sort -o "$OUT/$SAMPLE.trimmed.sorted.bam" "$OUT/$SAMPLE.trimmed.bam"
samtools index "$OUT/$SAMPLE.trimmed.sorted.bam"

# 3. Call variants: minimum depth 20, minimum base quality 20, no frequency
#    cut-off (-t 0). iVar's PASS column reflects its Fisher's exact test (p < 0.05).
samtools mpileup -aa -A -B -d 0 -Q 0 --reference "$REFERENCE" "$OUT/$SAMPLE.trimmed.sorted.bam" \
  | ivar variants -p "$OUT/ivar_tsv/$SAMPLE" -q 20 -t 0 -m 20 -r "$REFERENCE"

# 4. Per-position pileup for depth and A/C/G/T counts (Q20 is applied by the
#    Python converter), and depth for genome coverage at 10x.
samtools mpileup -aa -A -B -d 0 -Q 0 --reference "$REFERENCE" "$OUT/$SAMPLE.trimmed.sorted.bam" \
  > "$OUT/$SAMPLE.mpileup"
samtools depth -a "$OUT/$SAMPLE.trimmed.sorted.bam" > "$OUT/$SAMPLE.depth.tsv"

# 5. Convert to the website's input tables.
python offline_workflow/convert_offline_outputs.py pileup "$OUT/$SAMPLE.mpileup" "$SAMPLE" "$OUT/${SAMPLE}_base_counts.csv"
echo "Done. After all samples: python offline_workflow/convert_offline_outputs.py ivar $OUT/ivar_tsv $OUT/pooled_snvs.tsv"
