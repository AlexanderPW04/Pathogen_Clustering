"""Write the SYNTHETIC technical-replicate demonstration files (Sprint 5 Phase 2).

Everything here is invented: the 600-base reference "SYNREP.1", the sample and
specimen names, dates, site, read counts, and the true overdispersion values.
It is not real surveillance data and must not be cited as evidence about any
pathogen or sequencing protocol.

Read counts are drawn from a beta-binomial model with a known rho for each
specimen (seed 20261115), so the calibration result can be checked:

- SPEC_A: two replicates, true rho = 0.01, 40 variable sites plus one
  position with two alternative alleles (excluded from calibration).
- SPEC_B: two replicates, true rho = 0.002, 30 variable sites.
- SPEC_C: two replicates, but only 6 variable sites (too few to estimate).
- SPEC_D: two samples share a Specimen_ID but were collected a week apart,
  so they are NOT technical replicates and are not used.
- WW_SINGLE: one sample with no Specimen_ID (no replicate).

Every sample has its own Library_ID (LIB_<Sample_ID>), added in Phase 3: a
shared Specimen_ID counts as technical replication only when the two samples
document distinct libraries. The random draws are unchanged by that column.

A simulated variant caller reports an allele only when it has at least 2
reads and a frequency of at least 2%, so some low-frequency alleles are
reported in one replicate only; the A/C/G/T table (every position of every
sample) supplies the other replicate's count.

Run from the project folder:

    python examples/integrated/make_replicate_example.py
"""

from pathlib import Path

import numpy as np
import pandas as pd


FOLDER = Path(__file__).resolve().parent
SEED = 20261115
REFERENCE_ID = "SYNREP.1"
LENGTH = 600
LOCATION = "Synthetic Replicate Site"
BASES = "ACGT"
REPORT_MIN_READS = 2
REPORT_MIN_FREQUENCY = 0.02

# specimen, samples, dates, true rho, number of variable sites, adds a multiallelic site
SPECIMENS = [
    ("SPEC_A", ["REP_A1", "REP_A2"], ["2026-03-02", "2026-03-02"], 0.01, 40, True),
    ("SPEC_B", ["REP_B1", "REP_B2"], ["2026-03-09", "2026-03-09"], 0.002, 30, False),
    ("SPEC_C", ["REP_C1", "REP_C2"], ["2026-03-16", "2026-03-16"], 0.01, 6, False),
    ("SPEC_D", ["REP_D1", "REP_D2"], ["2026-03-23", "2026-03-30"], 0.01, 20, False),
    ("", ["WW_SINGLE"], ["2026-04-06"], 0.01, 25, False),
]
MULTIALLELIC_POSITION = 597


def reference_sequence(generator):
    return "".join(generator.choice(list(BASES), LENGTH))


def draw(generator, depth, p, rho):
    alpha, beta = p * (1 - rho) / rho, (1 - p) * (1 - rho) / rho
    return int(generator.binomial(depth, generator.beta(alpha, beta)))


def build(generator):
    reference = reference_sequence(generator)
    # A shared pool of variable positions (so mutations recur across dates) with fixed Alt bases.
    pool = sorted(generator.choice(np.arange(20, 580), 45, replace=False).tolist())
    alt_for = {position: generator.choice([b for b in BASES if b != reference[position - 1]])
               for position in pool}
    snv_rows, count_rows = [], []
    for specimen, samples, dates, rho, sites, multiallelic in SPECIMENS:
        positions = sorted(generator.choice(pool, sites, replace=False).tolist())
        # Mostly 3-60%, with some low-frequency alleles that the caller misses in one replicate.
        frequency = {position: float(generator.choice([generator.uniform(0.03, 0.6), generator.uniform(0.01, 0.03)],
                                                       p=[0.8, 0.2])) for position in positions}
        for sample, date in zip(samples, dates):
            library_depth = int(generator.integers(400, 2500))
            alt_reads = {}
            depth_at = {}
            for position in range(1, LENGTH + 1):
                depth = max(150, int(library_depth * generator.uniform(0.6, 1.4)))
                depth_at[position] = depth
                counts = dict.fromkeys(BASES, 0)
                ref = reference[position - 1]
                if position in frequency:
                    alt = draw(generator, depth, frequency[position], rho)
                    counts[alt_for[position]] = alt
                    counts[ref] = depth - alt
                    alt_reads[position] = [(alt_for[position], alt)]
                elif multiallelic and position == MULTIALLELIC_POSITION:
                    others = [b for b in BASES if b != ref][:2]
                    first, second = draw(generator, depth, 0.2, rho), draw(generator, depth, 0.1, rho)
                    second = min(second, depth - first)
                    counts[others[0]], counts[others[1]] = first, second
                    counts[ref] = depth - first - second
                    alt_reads[position] = [(others[0], first), (others[1], second)]
                else:
                    counts[ref] = depth
                count_rows.append({"Sample_ID": sample, "Reference_ID": REFERENCE_ID, "Position": position,
                                   "A_Count": counts["A"], "C_Count": counts["C"], "G_Count": counts["G"],
                                   "T_Count": counts["T"]})
            for position, alleles in sorted(alt_reads.items()):
                for alt, reads in alleles:
                    depth = depth_at[position]
                    if reads >= REPORT_MIN_READS and reads / depth >= REPORT_MIN_FREQUENCY:
                        snv_rows.append({"Sample_ID": sample, "Specimen_ID": specimen, "Library_ID": f"LIB_{sample}",
                                         "Reference_Genome": REFERENCE_ID, "Position": position,
                                         "Ref": reference[position - 1], "Alt": alt,
                                         "Alt_Frequency": round(reads / depth, 6), "Read_Depth": depth,
                                         "Alt_Depth": reads, "Collection_Date": date, "Location_Name": LOCATION})
    return reference, pd.DataFrame(snv_rows), pd.DataFrame(count_rows)


def write_all(folder=FOLDER):
    generator = np.random.default_rng(SEED)
    reference, snvs, counts = build(generator)
    folder = Path(folder)
    lines = [reference[start:start + 70] for start in range(0, LENGTH, 70)]
    (folder / "synthetic_replicate_reference.fasta").write_text(
        f">{REFERENCE_ID} SYNTHETIC invented reference for the technical-replicate demonstration "
        "(not a real organism)\n" + "\n".join(lines) + "\n", encoding="utf-8")
    snvs.to_csv(folder / "synthetic_replicate_snvs.csv", index=False)
    counts.to_csv(folder / "synthetic_replicate_base_counts.csv", index=False)
    return snvs, counts


if __name__ == "__main__":
    written, _ = write_all()
    print(f"Wrote {len(written)} SYNTHETIC replicate SNV records to {FOLDER}")
