"""Compare this application's import of Fontenele et al. (2021) Supplementary
Table 2 (mmc2.xlsx) with the counts the paper reports.

The workbook is not stored in this repository. Download it yourself from the
article's supplementary material (PMC8464352; the Europe PMC REST endpoint
https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8464352/supplementaryFiles
returns a ZIP that contains it), then run:

    python validation/check_fontenele_2021.py path/to/mmc2.xlsx

Every "paper" value below is quoted from the article text (Sections 3.2-3.3);
every "application" value is calculated here from the workbook. A difference
is reported, never hidden or adjusted.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomics.integrated_analysis import run_integrated_analysis  # noqa: E402
from wastewater.processed_snv import PAPER_FORMAT, load_processed_snvs  # noqa: E402


def compare(path):
    """Return a table of paper-reported versus application-calculated values."""

    data, warnings = load_processed_snvs(path)
    if data.attrs.get("input_format") != PAPER_FORMAT:
        raise ValueError("This file was not recognized as the published Supplementary Table 2.")
    distinct = data.drop_duplicates(["Position", "Ref", "Alt"])
    per_sample = data.groupby("Sample_ID").size()
    states = data["Location_Name"].str.split(":").str[0].str.strip()
    novel = data[data["Published_Novel_Global"].astype(str).str.strip().str.casefold().eq("novel")]
    novel_samples = novel.groupby(["Position", "Ref", "Alt"])["Sample_ID"].nunique()
    d614g = data[data["Position"].eq(23403) & data["Ref"].eq("A") & data["Alt"].eq("G")]
    low_depth = data["Read_Depth"].lt(20).sum()
    syn = distinct["Published_Syn_NonSyn"].astype(str).str.strip().value_counts()

    result = run_integrated_analysis(snv_path=path, input_labels={"snv": Path(path).name})
    summary = result["summary"]

    rows = [
        ("SNV records (all samples)", "7,973", len(data)),
        ("Samples with >90% coverage (analyzed)", "52", data["Sample_ID"].nunique()),
        ("States represented", "10", states.nunique()),
        ("SNVs per sample (range)", "24 to 793", f"{per_sample.min()} to {per_sample.max()}"),
        ("Distinct SNVs (position, Ref, Alt)", "5,680 unique SNVs", len(distinct)),
        ("Non-synonymous / synonymous distinct SNVs", "4,372 / 1,308",
         f"{syn.get('Nsyn', 0)} / {syn.get('Syn', 0)} (published 'Sy/Ns' column)"),
        ("Samples with spike D614G (A23403G)", "all 52", d614g["Sample_ID"].nunique()),
        ("Distinct SNVs flagged novel (any number of samples)", "not stated",
         len(novel_samples)),
        ("Novel SNVs seen in two or more samples", "548 (singletons removed)",
         int(novel_samples.ge(2).sum())),
        ("Novel SNVs in more than 8 samples", "4 (all non-synonymous)",
         int(novel_samples.gt(8).sum())),
        ("Collection dates of the 52 samples", "7 Apr to 16 Jun 2020 (all 91 samples)",
         f"{data['Collection_Date'].min().date()} to {data['Collection_Date'].max().date()}"),
        ("Records below the paper's own 20x depth minimum", "0 expected",
         int(low_depth)),
        ("Records retained by this app's paper-mode QC", "-", summary["pooled_records_included"]),
        ("Population structure available from this table", "Yue-Clayton PCoA (needs A/C/G/T counts)",
         summary.get("population_method")),
    ]
    table = pd.DataFrame(rows, columns=["Quantity", "Paper reports", "Application calculates"])
    return table, warnings


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("workbook", type=Path, help="Path to mmc2.xlsx")
    parser.add_argument("--output", type=Path, help="Optional CSV for the comparison table")
    arguments = parser.parse_args()
    table, warnings = compare(arguments.workbook)
    with pd.option_context("display.max_colwidth", 60, "display.width", 160):
        print(table.to_string(index=False))
    for warning in warnings:
        print("Import note:", warning)
    if arguments.output:
        table.to_csv(arguments.output, index=False)
        print("Saved", arguments.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
