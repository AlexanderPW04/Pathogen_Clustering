"""Write the SYNTHETIC mutation-mixture demonstration files (Sprint 5 Phase 3).

Everything here is invented: the 1,000-base reference "SYNMIX.1", the sites,
dates, specimens, libraries, read counts, and the "true" groups. It is not
real surveillance data and must not be cited as evidence about any pathogen.

The frequencies come from the Phase 3 simulation model
(validation/simulate_mutation_groups.py, planted structure), seed 20261115:

- 2 sites x 6 weekly dates = 12 specimens;
- 3 latent sources with logistic trajectories; groups of 10, 6, and 3
  mutations follow one source each;
- 2 "shared" mutations carry the sum of two sources, and 8 background
  mutations follow their own trajectories (none of these 10 belongs to a
  group);
- every specimen has two separately prepared libraries (distinct Library_IDs)
  with about 300 reads each and true overdispersion rho = 0.005, so the
  website can estimate rho and run the empirically calibrated mode;
- a caller reports an allele with at least 2 reads and frequency at least 2%;
  the A/C/G/T table covers every position of every library, so unreported
  mutations are measured non-detections, not unknowns.

`synthetic_mixture_truth.csv` lists the planted group of each mutation so the
result can be checked. It is not an input to the website.

Run from the project folder:

    python examples/integrated/make_mixture_example.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation.simulate_mutation_groups import beta_binomial, true_frequencies  # noqa: E402

FOLDER = Path(__file__).resolve().parent
SEED = 20261115
REFERENCE_ID = "SYNMIX.1"
LENGTH = 1000
SITES = ["Site North", "Site South"]
DATES = pd.date_range("2026-02-02", periods=6, freq="7D")
RHO = 0.005
LIBRARY_DEPTH = 300
BASES = "ACGT"
REPORT_MIN_READS, REPORT_MIN_FREQUENCY = 2, 0.02


def build(generator):
    reference = "".join(generator.choice(list(BASES), LENGTH))
    p, truth = true_frequencies(generator, "planted", sites=len(SITES), dates=len(DATES))
    count = len(truth)
    positions = np.sort(generator.choice(np.arange(20, LENGTH - 20), count, replace=False))
    alts = [generator.choice([b for b in BASES if b != reference[pos - 1]]) for pos in positions]
    units = [(site, date) for site in SITES for date in DATES]
    snv_rows, count_rows = [], []
    for unit_index, (site, date) in enumerate(units):
        specimen = f"MIX_{site.split()[1][0]}{unit_index % len(DATES) + 1}"
        for library in (1, 2):
            sample = f"{specimen}_L{library}"
            depth = np.maximum(150, np.round(LIBRARY_DEPTH * generator.lognormal(0, 0.25, LENGTH))).astype(int)
            alt_reads = {}
            for index, position in enumerate(positions):
                n = np.array([float(depth[position - 1])])
                alt_reads[position] = int(beta_binomial(generator, n, np.array([p[index, unit_index]]), RHO)[0])
            for position in range(1, LENGTH + 1):
                counts = dict.fromkeys(BASES, 0)
                ref = reference[position - 1]
                total = int(depth[position - 1])
                if position in alt_reads:
                    alt = alts[int(np.flatnonzero(positions == position)[0])]
                    counts[alt] = alt_reads[position]
                    counts[ref] = total - alt_reads[position]
                else:
                    counts[ref] = total
                count_rows.append({"Sample_ID": sample, "Reference_ID": REFERENCE_ID, "Position": position,
                                   "A_Count": counts["A"], "C_Count": counts["C"], "G_Count": counts["G"],
                                   "T_Count": counts["T"]})
            for index, position in enumerate(positions):
                reads, total = alt_reads[position], int(depth[position - 1])
                if reads >= REPORT_MIN_READS and reads / total >= REPORT_MIN_FREQUENCY:
                    snv_rows.append({"Sample_ID": sample, "Specimen_ID": specimen, "Library_ID": f"LIB_{sample}",
                                     "Reference_Genome": REFERENCE_ID, "Position": int(position),
                                     "Ref": reference[position - 1], "Alt": alts[index],
                                     "Alt_Frequency": round(reads / total, 6), "Read_Depth": total,
                                     "Alt_Depth": reads, "Collection_Date": date.date().isoformat(),
                                     "Location_Name": site})
    labels = {0: "Planted group A (10)", 1: "Planted group B (6)", 2: "Planted group C (3)"}
    truth_rows = [{"Variant_Key": f"{REFERENCE_ID}:{int(pos)}:{reference[pos - 1]}>{alt}",
                   "Planted_Group": labels.get(int(label), "None (shared)" if 100 <= label < 200 else "None (background)")}
                  for pos, alt, label in zip(positions, alts, truth)]
    return reference, pd.DataFrame(snv_rows), pd.DataFrame(count_rows), pd.DataFrame(truth_rows)


def write_all(folder=FOLDER):
    generator = np.random.default_rng(SEED)
    reference, snvs, counts, truth = build(generator)
    folder = Path(folder)
    lines = [reference[start:start + 70] for start in range(0, LENGTH, 70)]
    (folder / "synthetic_mixture_reference.fasta").write_text(
        f">{REFERENCE_ID} SYNTHETIC invented reference for the mutation-grouping demonstration "
        "(not a real organism)\n" + "\n".join(lines) + "\n", encoding="utf-8")
    snvs.to_csv(folder / "synthetic_mixture_snvs.csv", index=False)
    counts.to_csv(folder / "synthetic_mixture_base_counts.csv", index=False)
    truth.to_csv(folder / "synthetic_mixture_truth.csv", index=False)
    return snvs, counts, truth


if __name__ == "__main__":
    written, _, _ = write_all()
    print(f"Wrote {len(written)} SYNTHETIC mixture SNV records to {FOLDER}")
