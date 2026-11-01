# Optional offline raw-read workflow

This folder documents how raw paired-end reads could be turned into the
processed tables that the website analyzes. **It is optional and was not run
in this repository's environment**: BWA, SAMtools, and iVar are not installed
here, and raw FASTQ files are far too large for the Render Free web service.
The website never requires it; uploading processed tables is the supported
path.

| Stage | Tool (as in Fontenele et al. 2021) | Status in this project |
|---|---|---|
| Map reads to MN908947 / NC_045512.2 | `bwa mem` | Command template only (`raw_reads_to_tables.sh`) |
| Primer soft-clipping and quality trimming (window 4, PHRED 20, reads < 30 nt removed) | `ivar trim` | Command template only |
| Variant calling: depth >= 20, base quality >= 20, no frequency cut-off, Fisher p < 0.05 | `samtools mpileup` + `ivar variants` | Command template only |
| Per-position depth and coverage at 10x | `samtools depth -a`, `samtools mpileup -aa` | Command template only |
| Combine per-sample iVar TSVs (adds SAMPLE) | `convert_offline_outputs.py ivar` | Implemented and unit-tested on synthetic text |
| mpileup -> A/C/G/T counts with Q20 | `convert_offline_outputs.py pileup` | Implemented and unit-tested on synthetic text |
| Problematic-site masking (De Maio et al. 2020) | site list | Not bundled; upload a list on the website |

The primer BED file depends on the amplicon panel used for library
preparation and must come from that protocol. The paper's raw data are in
SRA BioProject PRJNA662596.

The ZooPhy/radx repository shows a similar BWA/SAMtools/iVar ordering, but it
has no license file, so no code was copied from it.

## Converting outputs

```bash
python offline_workflow/convert_offline_outputs.py ivar OUTPUT/ivar_tsv pooled_snvs.tsv
python offline_workflow/convert_offline_outputs.py pileup OUTPUT/S1.mpileup S1 S1_base_counts.csv
```

Concatenate the per-sample base-count CSVs, then upload `pooled_snvs.tsv` and
the combined counts (or run `run_integrated_analysis.py` locally when the
files exceed the 10 MB website limit).
