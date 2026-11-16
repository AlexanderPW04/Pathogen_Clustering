"""Phase 4A real-data validation of the generalized importer on published wastewater datasets.

Each dataset is imported through the same functions the website uses
(`inspect_file` -> `review_dataset` -> `prepare_inputs` ->
`run_integrated_analysis`), with the choices a user would confirm on the
review page. The imported records are then compared with values computed
*independently* from the original files (openpyxl / pandas only, no
application code), and with the published definitions.

The published files are only read, never changed; md5 checksums are printed
before and after.

Usage (each dataset is optional; missing folders are reported, not faked):

    python validation/check_published_datasets.py --output <folder>
        [--fontenele <mmc2.xlsx>]
        [--crits-christoph <folder with mBio.02703-20-st001.xlsx, -st004.xlsx,
                            reference.fasta, inStrain/*_SNVs.tsv>]
        [--schumann <folder with mmc1.xlsx, mmc2.xlsx, mmc3.xlsx>]

Retrieval (all public):
  * Fontenele et al. 2021, Water Research, doi:10.1016/j.watres.2021.117710:
    Supplementary Table 2 (mmc2.xlsx) from the article page.
  * Crits-Christoph et al. 2021, mBio, doi:10.1128/mBio.02703-20: Tables S1-S4 from
    https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7845645/supplementaryFiles ;
    reference.fasta and data/wastewater/<sample>/output/<sample>_SNVs.tsv from
    https://github.com/alexcritschristoph/wastewater_sarscov2 (GPL-3.0).
  * Schumann et al. 2022, Sci Total Environ, doi:10.1016/j.scitotenv.2022.158931:
    mmc1-mmc3 from
    https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9549760/supplementaryFiles .
"""

import argparse
import hashlib
import json
import re
import shutil
import sys
import tempfile
import time
from pathlib import Path

import openpyxl
import pandas as pd
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomics.input_inspection import inspect_file, prepare_inputs, review_dataset, suggest_roles  # noqa: E402
from genomics.integrated_analysis import run_integrated_analysis  # noqa: E402

NCBI_FASTA = ROOT / "examples" / "annotation" / "NC_045512.2.fasta"
NCBI_GENBANK = ROOT / "examples" / "annotation" / "NC_045512.2.gb"
SNV_HEADER = re.compile(r"::([ACGT])([1-9]\d*)([ACGT])$")


def md5(path):
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def run_import(files, options=None, settings=None):
    """files: [{"path", "role", "sheet"?, "options"?}] -> (review, result or None, seconds)."""

    start = time.perf_counter()
    work = Path(tempfile.mkdtemp())
    staged, infos = {}, []
    for index, item in enumerate(files):
        copy = work / str(index) / Path(item["path"]).name
        copy.parent.mkdir()
        shutil.copyfile(item["path"], copy)
        staged[index] = copy
        file_options = item.get("options", {})
        infos.append(inspect_file(copy, copy.name, item.get("sheet"), file_options.get("header_row"),
                                  file_options.get("sample_from_filename", False)))
    suggest_roles(infos)
    options = dict(options or {})
    options["files"] = {index: item.get("options", {}) for index, item in enumerate(files)}
    roles = {index: item["role"] for index, item in enumerate(files)}
    review = review_dataset(staged, infos, roles, {}, options)
    result = None
    if review["ready"]:
        folder = work / "prepared"
        folder.mkdir()
        paths, names = prepare_inputs(staged, review, folder, options)
        notes, provenance = names.pop("_notes", []), names.pop("_provenance", [])
        result = run_integrated_analysis(
            reference_path=paths.get("reference"), snv_path=paths.get("snv"),
            pooled_metadata_path=paths.get("pooled_metadata"), position_path=paths.get("positions"),
            annotation_path=paths.get("annotation"),
            settings={"qc_mode": "paper", **(settings or {})}, input_labels=names, input_notes=notes,
            import_details={"provenance": provenance, "readiness": review["analysis_readiness"],
                            "matching": review.get("matching")})
    shutil.rmtree(work, ignore_errors=True)
    return review, result, round(time.perf_counter() - start, 1)


def check(results, name, passed, detail):
    results.append({"Check": name, "Result": "PASS" if passed else "FAIL", "Detail": detail})


def readiness_summary(review):
    return {row["Analysis"]: row["Status"] for row in review["analysis_readiness"]}


def reference_sequence(path):
    return str(next(SeqIO.parse(str(path), "fasta")).seq).upper()


# ---------------------------------------------------------------------------
# Dataset A: Fontenele et al. (2021) mmc2.xlsx (regression dataset)
# ---------------------------------------------------------------------------

def fontenele(path, output):
    results = []
    before = md5(path)
    review, result, seconds = run_import([{"path": path, "role": "snv"}, {"path": NCBI_FASTA, "role": "reference"}])
    table = result["tables"]["canonical_variant_observations.csv"]
    samples = result["tables"]["canonical_samples.csv"]
    workbook = pd.read_excel(path, header=1, dtype=str)  # independent read: title row, then headers
    check(results, "All published records imported", len(table) == len(workbook),
          f"{len(table)} imported vs {len(workbook)} rows in the workbook")
    check(results, "Sample count", samples["Sample_ID"].nunique() == workbook["SAMPLE"].nunique(),
          f"{samples['Sample_ID'].nunique()} samples vs {workbook['SAMPLE'].nunique()}")
    reported = pd.to_numeric(table["Reported_VAF"]).to_numpy()
    published = pd.to_numeric(workbook["ALT_FREQ"]).to_numpy()
    check(results, "Reported frequencies unchanged (fractions)", (abs(reported - published) < 1e-12).all(),
          f"max difference {abs(reported - published).max():.2e}; published max {published.max()}")
    depth = pd.to_numeric(table["Total_Depth"]).to_numpy()
    check(results, "Depths equal TOTAL_DP", (depth == pd.to_numeric(workbook["TOTAL_DP"]).to_numpy()).all(),
          "TOTAL_DP is the total depth (iVar definition)")
    refs = result["tables"]["reference_verification.csv"]
    check(results, "Reference bases verified against the NCBI NC_045512.2 FASTA",
          refs["Status"].eq("Verified compatible").all() and int(refs["Ref_Base_Mismatches"].fillna(0).sum()) == 0,
          refs[["Reference_ID", "Status", "Ref_Base_Matches", "Ref_Base_Mismatches"]].to_dict("records"))
    measured = table["Count_Provenance"].value_counts().to_dict()
    check(results, "Alternative counts are reconstructed, never labelled measured",
          "Supplied" not in measured, measured)
    after = md5(path)
    check(results, "Original file unchanged", before == after, before)
    return {"dataset": "Fontenele et al. 2021 (mmc2.xlsx)", "seconds": seconds, "checks": results,
            "readiness": readiness_summary(review), "records": len(table),
            "samples": int(samples["Sample_ID"].nunique())}


# ---------------------------------------------------------------------------
# Dataset B: Crits-Christoph et al. (2021)
# ---------------------------------------------------------------------------

def expected_alleles(row, minimum):
    return [base for base in "ACGT" if base != row["ref_base"] and int(row[base]) >= minimum]


def crits_christoph(folder, output):
    folder = Path(folder)

    def find(name):
        found = sorted(folder.rglob(name))
        if not found:
            raise FileNotFoundError(f"{name} was not found under {folder}")
        return found[0]

    s1, s4 = find("mBio.02703-20-st001.xlsx"), find("mBio.02703-20-st004.xlsx")
    reference = find("reference.fasta")
    checksums = {path.name: md5(path) for path in (s1, s4, reference)}
    results, runs = [], {}
    label = next(SeqIO.parse(str(reference), "fasta")).id
    s4_options = {"base_count_minimum": 2, "declared_reference": label}
    # 1. Table S4 with Table S1 metadata: the join must be refused (IDs differ), then accepted explicitly.
    review, _, _ = run_import([{"path": s4, "role": "snv", "options": s4_options},
                               {"path": reference, "role": "reference"},
                               {"path": s1, "role": "pooled_metadata"}])
    matching = {row["Check"]: row for row in review["matching"]}
    matched = matching["Variant samples matched to a metadata row"]["Count"]
    check(results, "S4-S1 join refused without identical IDs",
          not review["ready"] and matched == 0,
          f"{matched} of {len(review['pooled_ids'])} S4 samples have an S1 row; problem shown: "
          + next((text for text in review["problems"] if "missing from the pooled metadata" in text), "none"))
    check(results, "S1 summary row without Sample_ID reported",
          "Metadata rows without a Sample_ID" in matching, matching.get("Metadata rows without a Sample_ID"))
    review, result, seconds = run_import([{"path": s4, "role": "snv", "options": s4_options},
                                          {"path": reference, "role": "reference"},
                                          {"path": s1, "role": "pooled_metadata"}],
                                         options={"allow_unmatched_metadata": True})
    runs["S4 + S1 (unmatched accepted)"] = {"seconds": seconds, "readiness": readiness_summary(review)}
    observations = result["tables"]["canonical_variant_observations.csv"]
    samples = result["tables"]["canonical_samples.csv"]
    check(results, "No dates were joined from S1", samples["Collection_Date"].isna().all()
          or samples["Collection_Date"].astype(str).isin(["", "nan", "NaT"]).all(),
          "Dates stay unknown because no ID matches")
    # Independent derivation from the raw table.
    raw = pd.read_excel(s4, dtype=str)
    expected = []
    for _, row in raw.iterrows():
        depth = sum(int(row[base]) for base in "ACGT")
        for base in expected_alleles(row, 2):
            expected.append((row["sample"], int(row["genome_position"]), row["ref_base"], base,
                             int(row[base]), depth))
    imported = list(zip(observations["Sample_ID"], observations["Genome_Position"].astype(int),
                        observations["Ref"], observations["Alt"], observations["Alt_Read_Count"].astype(int),
                        observations["Total_Depth"].astype(int)))
    check(results, "Records equal an independent derivation (>= 2 reads, the paper's stated minimum)",
          sorted(expected) == sorted(imported), f"{len(imported)} imported, {len(expected)} expected")
    # Agreement with the study's own reported allele and frequency (con_base / var_base, con_freq / var_freq).
    agree = differ = 0
    extra = []
    key = observations.set_index(["Sample_ID", "Genome_Position", "Alt"])["Calculated_VAF"]
    for _, row in raw.iterrows():
        allele = row["con_base"] if row["con_base"] != row["ref_base"] else row["var_base"]
        frequency = float(row["con_freq"] if row["con_base"] != row["ref_base"] else row["var_freq"])
        index = (row["sample"], int(row["genome_position"]), allele)
        if allele != row["ref_base"] and index in key.index:
            agree += abs(float(key.loc[index]) - frequency) < 1e-9
            differ += abs(float(key.loc[index]) - frequency) >= 1e-9
        others = [base for base in expected_alleles(row, 2) if base != allele]
        extra += [f"{row['sample']}:{row['genome_position']}{row['ref_base']}>{base}" for base in others]
    check(results, "Calculated VAF equals the study's reported allele frequency", differ == 0 and agree > 0,
          f"{agree} of {len(raw)} rows agree to 1e-9; {differ} differ")
    check(results, "Extra alleles beyond the study's two reported bases are listed", True,
          f"{len(extra)} additional non-reference allele(s) with >= 2 reads: {', '.join(extra[:5])}")
    sequence = reference_sequence(reference)
    matches = sum(sequence[int(p) - 1] == r for p, r in zip(raw["genome_position"], raw["ref_base"]))
    check(results, "Positions are 1-based against WIV04 (independent)", matches == len(raw),
          f"{matches} of {len(raw)} Ref bases match at position - 1")
    refs = result["tables"]["reference_verification.csv"]
    check(results, "Application reference verification", refs["Status"].eq("Verified compatible").all(),
          refs[["Reference_ID", "Status", "Ref_Base_Matches"]].to_dict("records"))
    nc = reference_sequence(NCBI_FASTA)
    check(results, "WIV04 is a different sequence from NC_045512.2 (never treated as the same)",
          sequence != nc, f"WIV04 {len(sequence)} bases, NC_045512.2 {len(nc)}; prefix identical: "
          f"{nc.startswith(sequence)}")
    # 2. The GitHub inStrain per-sample files: 0-based positions and file-name sample IDs.
    instrain = sorted(folder.rglob("*_SNVs.tsv"))
    if instrain:
        files = [{"path": path, "role": "snv",
                  "options": {"sample_from_filename": True, "base_count_minimum": 2}} for path in instrain]
        files.append({"path": reference, "role": "reference"})
        review, _, _ = run_import(files)
        warned = [text for text in review["warnings"] if "match with 0-based positions" in text]
        check(results, "0-based numbering detected (warning, no automatic shift)", len(warned) == len(instrain),
              f"{len(warned)} of {len(instrain)} files warned")
        for item in files[:-1]:
            item["options"]["position_base"] = "0"
        review, result, seconds = run_import(files)
        runs["inStrain per-sample files (0-based confirmed)"] = {"seconds": seconds,
                                                                "readiness": readiness_summary(review)}
        converted = result["tables"]["canonical_variant_observations.csv"]
        converted = converted.assign(Sample=converted["Sample_ID"].str.replace("_SNVs", "", regex=False))
        shared = set(raw["sample"])
        sites = {(sample, int(position)) for sample, position in zip(converted["Sample"], converted["Genome_Position"])
                 if sample in shared}
        published = {(sample, int(position)) for sample, position in zip(raw["sample"], raw["genome_position"])}
        exact = len(sites & published)
        shifted = max(len({(a, b + shift) for a, b in sites} & published) for shift in (-1, 1))
        check(results, "inStrain per-sample files: converted positions coincide with Table S4 sites",
              exact > 5 * max(shifted, 1), f"{exact} of {len(sites)} sites at the same 1-based position (a one-base "
              f"error would give {shifted}); these files are an earlier inStrain run than S4: {len(sites - published)} "
              f"site(s) only in the files, {len(published - sites)} only in S4 (validation maps '<sample>_SNVs' "
              "file names to S4 sample names)")
    final = sorted(folder.rglob("Wastewater_SNVs_v3.csv"))
    if final:
        # The repository's final long table: 0-based positions, sample names with a pipeline suffix.
        files = [{"path": final[0], "role": "snv",
                  "options": {"base_count_minimum": 2, "position_base": "0"}},
                 {"path": reference, "role": "reference"}]
        review, result, seconds = run_import(files)
        runs["GitHub final_tables/Wastewater_SNVs_v3.csv (long, 0-based)"] = {
            "seconds": seconds, "readiness": readiness_summary(review), "problems": review["problems"]}
        if result is not None:
            table = result["tables"]["canonical_variant_observations.csv"]
            left = sorted(zip(table["Sample_ID"].str.split("-").str[0], table["Genome_Position"].astype(int),
                              table["Ref"], table["Alt"], table["Alt_Read_Count"].astype(int),
                              table["Total_Depth"].astype(int)))
            check(results, "0-based long table (+1) reproduces the Table S4 records exactly", left == sorted(expected),
                  f"{len(left)} vs {len(expected)} records (sample names compared before the '-mapq20...' suffix)")
    after = {path.name: md5(path) for path in (s1, s4, reference)}
    check(results, "Original files unchanged", after == checksums, checksums)
    return {"dataset": "Crits-Christoph et al. 2021 (Tables S1, S4; GitHub inStrain tables)", "checks": results,
            "runs": runs, "records": len(observations), "samples": int(samples["Sample_ID"].nunique())}


# ---------------------------------------------------------------------------
# Dataset C: Schumann et al. (2022)
# ---------------------------------------------------------------------------

def wide_cells(path, sheet, header_row):
    """Independent read of a wide sheet: {(sample, header): value text} for non-NA cells."""

    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    rows = workbook[sheet].iter_rows(values_only=True)
    for _ in range(header_row - 1):
        next(rows)
    header = [str(cell).strip() if cell is not None else "" for cell in next(rows)]
    cells, samples = {}, []
    for row in rows:
        if not row or row[0] is None:
            continue
        samples.append(str(row[0]))
        for name, value in zip(header, row):
            if SNV_HEADER.search(name) and value not in (None, "NA", "") and str(value).strip() != "NA":
                cells[(str(row[0]), name)] = value
    workbook.close()
    return cells, samples, header


def schumann(folder, output, include_large=True):
    folder = Path(folder)
    mmc1, mmc2, mmc3 = folder / "mmc1.xlsx", folder / "mmc2.xlsx", folder / "mmc3.xlsx"
    checksums = {path.name: md5(path) for path in (mmc1, mmc2, mmc3) if path.exists()}
    results, runs = [], {}
    sheet = "S2.1 - mutations_found_dataset-"
    files = [{"path": mmc1, "role": "snv", "sheet": sheet,
              "options": {"layout_confirmed": True, "declared_reference": "NC_045512.2"}},
             {"path": NCBI_FASTA, "role": "reference"}, {"path": NCBI_GENBANK, "role": "annotation"}]
    review, _, _ = run_import([dict(files[0], options={"declared_reference": "NC_045512.2"})])
    check(results, "Wide layout needs the user's confirmation", not review["ready"] and any(
        "Confirm the wide layout" in text for text in review["problems"]), "problem shown before confirmation")
    review, result, seconds = run_import(files)
    runs["S2.1 Berlin (38 samples, wide)"] = {"seconds": seconds, "readiness": readiness_summary(review)}
    observations = result["tables"]["canonical_variant_observations.csv"]
    cells, samples, header = wide_cells(mmc1, sheet, 2)
    check(results, "Every reported single-nucleotide cell became one record", len(observations) == len(cells),
          f"{len(observations)} records vs {len(cells)} non-NA single-nucleotide cells")
    lookup = {}
    for (sample, name), value in cells.items():
        ref, position, alt = SNV_HEADER.search(name).groups()
        lookup[(sample, int(position), ref, alt)] = float(value)
    exact = sum(abs(lookup.get((s, int(p), r, a), -1) - float(f)) < 1e-12 for s, p, r, a, f in zip(
        observations["Sample_ID"], observations["Genome_Position"], observations["Ref"], observations["Alt"],
        observations["Reported_VAF"]))
    check(results, "Each frequency equals its cell (sample and mutation assigned correctly)",
          exact == len(observations), f"{exact} of {len(observations)} exact")
    check(results, "NA cells became no record (never 0%)",
          not (pd.to_numeric(observations["Reported_VAF"]) == 0).any(),
          f"{len(samples) * sum(1 for name in header if SNV_HEADER.search(name)) - len(cells)} NA cells")
    check(results, "Depth is not invented", observations["Total_Depth"].isna().all(),
          observations["Measurement_Class"].value_counts().to_dict())
    counts = pd.read_excel(mmc1, sheet_name="S4 - mutation_counts", dtype=str)
    counts = counts[counts["sample"].ne("Total")].set_index("sample")["total_muts"].astype(float).astype(int)
    per_sample = observations.groupby("Sample_ID").size()
    other = {}
    workbook = openpyxl.load_workbook(mmc1, read_only=True, data_only=True)
    rows = list(workbook[sheet].iter_rows(values_only=True))
    workbook.close()
    names = [str(cell).strip() if cell is not None else "" for cell in rows[1]]
    first_mutation = names.index("coordinates_long") + 1
    for row in rows[2:]:
        if row and row[0]:
            other[str(row[0])] = sum(1 for name, value in zip(names[first_mutation:], row[first_mutation:])
                                     if not SNV_HEADER.search(name) and value not in (None, "NA"))
    agree = [sample for sample in counts.index
             if per_sample.get(sample, 0) + other.get(sample, 0) == counts[sample]]
    check(results, "Per-sample mutation counts match the authors' hidden 'S4 - mutation_counts' sheet",
          len(agree) == len(counts), f"{len(agree)} of {len(counts)} samples: imported SNV records + "
                                     "non-SNV columns = published total_muts")
    no_records = sorted(set(samples) - set(observations["Sample_ID"]))
    check(results, "Samples without any reported SNV are reported, not hidden",
          all(per_sample.get(sample, 0) == 0 for sample in no_records), ", ".join(no_records) or "none")
    dates = result["tables"]["canonical_samples.csv"].set_index("Sample_ID")["Collection_Date"].astype(str)
    raw_dates = pd.read_excel(mmc1, sheet_name=sheet, header=1, dtype=str).set_index("samplename")["dates"]
    same = sum(dates.get(sample, "")[:10] == str(value)[:10] for sample, value in raw_dates.items()
               if sample in dates.index)
    check(results, "Collection dates joined to the right samples", same == len(dates),
          f"{same} of {len(dates)} samples")
    refs = result["tables"]["reference_verification.csv"]
    check(results, "Ref bases agree with NC_045512.2 (PiGx default reference)",
          refs["Status"].eq("Verified compatible").all(), refs[["Reference_ID", "Status"]].to_dict("records"))
    # Phase 4 annotation versus the authors' VEP protein labels in the headers.
    annotated = result["tables"].get("sprint5_annotated_mutations.csv")
    if annotated is not None:
        ours = annotated.set_index("Nucleotide_Change")
        compared = agreed = 0
        disagreements = []
        for name in header:
            match = SNV_HEADER.search(name)
            if not match:
                continue
            ref, position, alt = match.groups()
            change = f"{position}{ref}>{alt}"
            protein = name.split("::")[0]
            gene, _, aa = protein.partition(":")
            vep = re.match(r"^([A-Z*])(\d+)([A-Z*-])$", aa)
            if change not in ours.index or not vep:
                continue
            compared += 1
            consequence = ours.at[change, "Amino_Acid_Changes"]
            expected = f"{vep.group(1)}{vep.group(2)}{vep.group(1) if vep.group(3) == '-' else vep.group(3)}"
            if isinstance(consequence, str) and expected in consequence:
                agreed += 1
            elif len(disagreements) < 8:
                disagreements.append(f"{name} vs {consequence}")
        check(results, "Phase 4 amino-acid predictions agree with the authors' VEP labels",
              compared > 0 and agreed / compared > 0.95,
              f"{agreed} of {compared} labelled SNVs agree; examples of differences: {disagreements}")
    if include_large and mmc2.exists():
        review, result, seconds = run_import([{"path": mmc2, "role": "snv",
                                               "options": {"layout_confirmed": True,
                                                           "declared_reference": "NC_045512.2"}},
                                              {"path": NCBI_FASTA, "role": "reference"}])
        if result is not None:
            runs["mmc2 Berlin35 (109 samples, 8,862 mutation columns)"] = {
                "seconds": seconds, "readiness": readiness_summary(review),
                "records": len(result["tables"]["canonical_variant_observations.csv"])}
            cells, _, _ = wide_cells(mmc2, "data_mutation_plot.csv", 1)
            check(results, "Berlin35: every reported cell became one record",
                  len(result["tables"]["canonical_variant_observations.csv"]) == len(cells),
                  f"{len(result['tables']['canonical_variant_observations.csv'])} vs {len(cells)}")
        else:
            runs["mmc2 Berlin35"] = {"problems": review["problems"]}
    if mmc3.exists():
        info = inspect_file(mmc3, mmc3.name)
        check(results, "UCSD mmc3 (14 MB) is refused with an actionable message, not a crash",
              info["refused"], "; ".join(info["problems"])[:300])
    after = {path.name: md5(path) for path in (mmc1, mmc2, mmc3) if path.exists()}
    check(results, "Original files unchanged", after == checksums, checksums)
    return {"dataset": "Schumann et al. 2022 (mmc1 S2.1, mmc2, mmc3)", "checks": results, "runs": runs,
            "records": len(observations)}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--output", required=True)
    parser.add_argument("--fontenele")
    parser.add_argument("--crits-christoph")
    parser.add_argument("--schumann")
    parser.add_argument("--skip-large", action="store_true", help="Skip the 109-sample Berlin35 matrix")
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    report = []
    for name, path, function in [("Fontenele", args.fontenele, fontenele),
                                 ("Crits-Christoph", args.crits_christoph, crits_christoph),
                                 ("Schumann", args.schumann, schumann)]:
        if not path or not Path(path).exists():
            report.append({"dataset": name, "checks": [], "note": "Not supplied: real-data checks not run."})
            continue
        if function is schumann:
            report.append(function(path, output, include_large=not args.skip_large))
        else:
            report.append(function(path, output))
    (output / "published_dataset_validation.json").write_text(json.dumps(report, indent=1, default=str))
    failed = 0
    for item in report:
        print("==", item["dataset"], item.get("note", ""))
        for row in item["checks"]:
            failed += row["Result"] == "FAIL"
            print(f"  {row['Result']}  {row['Check']}: {str(row['Detail'])[:240]}")
        for run, details in item.get("runs", {}).items():
            print(f"  run {run}: {details}")
    print("FAILED checks:", failed)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
