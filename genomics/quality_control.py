"""Quality control for pooled-sample variants, based on Fontenele et al. (2021).

The paper's methods (Sections 2.4-2.5) used:

- minimum base quality 20 and iVar's sliding-window quality trimming;
- genome coverage counted at 10x depth, with detailed analysis restricted to
  samples covering more than 90% of the genome;
- iVar variant calls with at least 20x depth, Fisher's exact test p < 0.05,
  and no minimum allele-frequency cut-off;
- removal of problematic sites suggested by De Maio et al. (2020).

A processed SNV table usually cannot prove every one of these steps. Each
filter is applied only when the needed values were supplied; otherwise it is
reported as "not evaluable" instead of being claimed as passed.

Two modes are available:

- "paper": records or samples that fail an evaluable paper filter are
  excluded from downstream analyses (they stay in the downloadable table).
- "exploratory": nothing is excluded; failures are only flagged.
"""

import pandas as pd

from genomics.pooled_observations import read_any_table, standard_name
from genomics.references import reference_base, references_compatible


PAPER_SETTINGS = {
    "min_base_quality": 20,
    "coverage_depth": 10,
    "min_genome_coverage_percent": 90.0,
    "min_variant_depth": 20,
    "max_fisher_p": 0.05,
    "min_alt_frequency": None,
}

PAPER_MODE = "paper"
EXPLORATORY_MODE = "exploratory"

PASS = "Pass"
NOT_EVALUABLE = "Not evaluable"


def load_masked_sites(path):
    """Read a list of problematic sites with Reference_ID and Position columns."""

    data = read_any_table(path)
    data.columns = [
        standard_name(name, {"reference_id": "Reference_ID", "position": "Position",
                             "pos": "Position", "chrom": "Reference_ID"})
        for name in data.columns
    ]
    if not {"Reference_ID", "Position"}.issubset(data.columns):
        raise ValueError("A masked-site list needs Reference_ID and Position columns.")
    positions = pd.to_numeric(data["Position"], errors="coerce")
    if positions.isna().any() or positions.lt(1).any() or positions.ne(positions.round()).any():
        raise ValueError("Masked-site positions must be whole numbers of 1 or more.")
    sites = set()
    for reference_id, position in zip(data["Reference_ID"], positions):
        sites.add((str(reference_id).strip(), int(position)))
    return sites


def depth_check(depth, minimum):
    if pd.isna(depth):
        return f"{NOT_EVALUABLE} (depth not supplied)"
    if depth >= minimum:
        return PASS
    return f"Fail (depth {int(depth)} < {minimum})"


def fisher_check(row, maximum):
    p_value = pd.to_numeric(row.get("Source_Fisher_P", pd.NA), errors="coerce")
    if pd.notna(p_value):
        if p_value < maximum:
            return PASS
        return f"Fail (p = {p_value:.3g}, not < {maximum})"
    caller_pass = str(row.get("Source_Caller_Pass", "")).strip().casefold()
    if caller_pass in {"true", "false"}:
        # iVar's PASS column is TRUE when its Fisher's exact test p < 0.05.
        return PASS if caller_pass == "true" else "Fail (caller PASS = FALSE)"
    return f"{NOT_EVALUABLE} (no Fisher p-value supplied)"


def coverage_check(percent, minimum):
    if pd.isna(percent):
        return f"{NOT_EVALUABLE} (genome coverage unknown)"
    if percent > minimum:
        return PASS
    return f"Fail ({percent:.1f}% genome coverage, not > {minimum:.0f}%)"


def reference_check(row, reference, registry):
    """Data-integrity check (not a paper step): does Ref match the reference?"""

    if reference is None:
        return "Not checked (no reference supplied)"
    ok, _ = references_compatible(row["Reference_ID"], reference["Reference_ID"], registry)
    if not ok:
        return "Not checked (record uses a different reference)"
    position = int(row["Position"])
    if position > reference["Length"]:
        return "Fail (position beyond the reference)"
    if reference_base(reference, position) != row["Ref"]:
        return "Fail (Ref base does not match the reference)"
    return PASS


def apply_pooled_quality_control(variants, samples, mode=PAPER_MODE, masked_sites=None,
                                 settings=None, reference=None, registry=None):
    """Return (variants, samples, method_report) with quality columns added."""

    if mode not in {PAPER_MODE, EXPLORATORY_MODE}:
        raise ValueError("Quality-control mode must be 'paper' or 'exploratory'.")
    settings = settings or PAPER_SETTINGS
    variants = variants.copy()
    samples = samples.copy()

    sample_results = []
    for row in samples.to_dict("records"):
        sample_results.append(
            coverage_check(row["Genome_Coverage_Percent"], settings["min_genome_coverage_percent"])
        )
    samples["QC_Genome_Coverage"] = sample_results
    coverage_by_sample = dict(zip(samples["Sample_ID"], samples["QC_Genome_Coverage"]))

    depth_results = []
    fisher_results = []
    mask_results = []
    reference_results = []
    for row in variants.to_dict("records"):
        reference_results.append(reference_check(row, reference, registry))
        depth_results.append(depth_check(row["Read_Depth"], settings["min_variant_depth"]))
        fisher_results.append(fisher_check(row, settings["max_fisher_p"]))
        if masked_sites is None:
            mask_results.append("Not reproduced (no masked-site list supplied)")
        elif (row["Reference_ID"], int(row["Position"])) in masked_sites:
            mask_results.append("Fail (listed problematic site)")
        else:
            mask_results.append(PASS)

    variants["QC_Variant_Depth"] = depth_results
    variants["QC_Fisher_Test"] = fisher_results
    variants["QC_Masked_Site"] = mask_results
    variants["QC_Base_Quality"] = "Not verifiable from processed records"
    variants["QC_Min_Frequency"] = "No cut-off (as in the paper)"
    variants["QC_Sample_Coverage"] = variants["Sample_ID"].map(coverage_by_sample)
    variants["QC_Reference_Base"] = reference_results

    check_columns = ["QC_Variant_Depth", "QC_Fisher_Test", "QC_Masked_Site",
                     "QC_Sample_Coverage", "QC_Reference_Base"]
    statuses = []
    included = []
    for row in variants[check_columns].to_dict("records"):
        failures = [name.replace("QC_", "").replace("_", " ").lower()
                    for name, value in row.items() if str(value).startswith("Fail")]
        if not failures:
            statuses.append("Retained")
            included.append(True)
        elif mode == PAPER_MODE:
            statuses.append("Excluded by paper preset: " + ", ".join(failures))
            included.append(False)
        else:
            statuses.append("Flagged (exploratory, retained): " + ", ".join(failures))
            included.append(True)
    variants["Quality_Status"] = statuses
    variants["Analysis_Included"] = included

    sample_status = []
    for value in samples["QC_Genome_Coverage"]:
        if str(value).startswith("Fail"):
            if mode == PAPER_MODE:
                sample_status.append("Excluded by paper preset: genome coverage")
            else:
                sample_status.append("Flagged (exploratory, retained): genome coverage")
        else:
            sample_status.append("Retained")
    samples["Quality_Status"] = sample_status
    samples["Analysis_Included"] = [not status.startswith("Excluded") for status in sample_status]

    report = method_report(variants, samples, masked_sites, settings, mode)
    return variants, samples, report


def evaluable_count(values):
    return int(sum(not str(value).startswith(NOT_EVALUABLE) for value in values))


def method_report(variants, samples, masked_sites, settings, mode):
    """Describe which paper methods were reproduced, applied, or unavailable."""

    total_variants = len(variants)
    total_samples = len(samples)
    depth_evaluable = evaluable_count(variants["QC_Variant_Depth"])
    fisher_evaluable = evaluable_count(variants["QC_Fisher_Test"])
    coverage_evaluable = evaluable_count(samples["QC_Genome_Coverage"])

    def applied_status(evaluable, total):
        if total and evaluable == total:
            return "Applied to supplied values"
        if evaluable:
            return "Partly applied (some values missing)"
        return "Not evaluable"

    rows = [
        {"Step": "Read mapping to the reference",
         "Paper_Setting": "BWA-MEM to MN908947 (RefSeq NC_045512.2)",
         "Status": "Not reproduced here",
         "Detail": "Processed records were uploaded; reads were not re-mapped. See the optional offline workflow."},
        {"Step": "Primer and quality trimming",
         "Paper_Setting": "iVar trim; sliding window 4, minimum PHRED 20; reads < 30 nt removed",
         "Status": "Not reproduced here",
         "Detail": "Requires raw reads and the primer scheme."},
        {"Step": "Minimum base quality",
         "Paper_Setting": f"Q{settings['min_base_quality']}",
         "Status": "Not verifiable",
         "Detail": "Processed SNV records do not show the base qualities that were counted."},
        {"Step": "Sample genome coverage",
         "Paper_Setting": (f">{settings['min_genome_coverage_percent']:.0f}% of the genome "
                           f"at >= {settings['coverage_depth']}x"),
         "Status": applied_status(coverage_evaluable, total_samples),
         "Detail": f"Genome coverage known for {coverage_evaluable} of {total_samples} pooled samples."},
        {"Step": "Minimum variant depth",
         "Paper_Setting": f">= {settings['min_variant_depth']}x",
         "Status": applied_status(depth_evaluable, total_variants),
         "Detail": f"Depth supplied for {depth_evaluable} of {total_variants} variant records."},
        {"Step": "Fisher's exact test",
         "Paper_Setting": f"iVar p < {settings['max_fisher_p']}",
         "Status": applied_status(fisher_evaluable, total_variants),
         "Detail": (f"p-values or caller PASS flags supplied for {fisher_evaluable} of "
                    f"{total_variants} records; p-values are never recalculated or invented."),
         },
        {"Step": "Minimum allele frequency",
         "Paper_Setting": "None (population-level data)",
         "Status": "Matched",
         "Detail": "No frequency cut-off is applied."},
        {"Step": "Problematic-site masking",
         "Paper_Setting": "De Maio et al. (2020) sites as of September 2020",
         "Status": "Applied to supplied list" if masked_sites is not None else "Not reproduced",
         "Detail": ("Positions in the uploaded list were checked." if masked_sites is not None
                    else "The historical site list is not bundled and was not supplied."),
         },
        {"Step": "Reference base agreement",
         "Paper_Setting": "Not a paper step (data-integrity check added here)",
         "Status": ("Checked against the supplied reference"
                    if not variants["QC_Reference_Base"].str.startswith("Not checked").all()
                    else "Not checked"),
         "Detail": "A Ref base that disagrees with the reference suggests a different coordinate system."},
        {"Step": "Quality-control mode",
         "Paper_Setting": "-",
         "Status": "Paper reproduction preset" if mode == PAPER_MODE else "Exploratory",
         "Detail": ("Records failing an evaluable filter are excluded from analyses."
                    if mode == PAPER_MODE else "Failures are flagged but retained."),
         },
    ]
    return pd.DataFrame(rows, columns=["Step", "Paper_Setting", "Status", "Detail"])
