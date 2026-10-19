import tempfile
import unittest
from pathlib import Path

import pandas as pd

from validate_inputs import load_and_validate_inputs


class InputValidationTests(unittest.TestCase):
    def valid_metadata(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "Accession": ["ACC001", "ACC002"],
                "Pathogen": ["Pathogen A", "Pathogen B"],
                "Gene_Segment": ["Segment 1", "Whole genome"],
                "Collection_Date": ["2026-01-15", "2026-02-20"],
                "Location_Name": ["Location A", "Location B"],
                "Latitude": [33.45, "40.71"],
                "Longitude": [-112.07, "-74.00"],
                "Host_Source": ["Human", "Environmental"],
            }
        )

    def write_fasta(self, folder: str, text: str) -> Path:
        path = Path(folder) / "sequences.fasta"
        path.write_text(text, encoding="utf-8")
        return path

    def write_metadata(
        self,
        folder: str,
        data: pd.DataFrame,
        extension: str,
    ) -> Path:
        path = Path(folder) / f"metadata{extension}"
        if extension == ".xlsx":
            data.to_excel(path, index=False, engine="openpyxl")
        else:
            data.to_csv(path, index=False)
        return path

    def test_valid_fasta_and_xlsx_are_matched_and_cleaned(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            fasta_path = self.write_fasta(
                folder,
                ">ACC001 description\nacgt\nnn--\n>ACC002\nACGT\n",
            )
            metadata_path = self.write_metadata(
                folder, self.valid_metadata(), ".xlsx"
            )
            data, errors, missing_counts = load_and_validate_inputs(
                fasta_path, metadata_path
            )

        self.assertEqual(errors, [])
        self.assertEqual(data["Sample_ID"].tolist(), ["ACC001", "ACC002"])
        self.assertEqual(data.loc[0, "Sequence"], "ACGTNN--")
        self.assertTrue(pd.api.types.is_datetime64_any_dtype(data["Collection_Date"]))
        self.assertTrue(pd.api.types.is_numeric_dtype(data["Latitude"]))
        self.assertTrue(all(count == 0 for count in missing_counts.values()))

    def test_csv_metadata_preserves_accession_text(self) -> None:
        metadata = self.valid_metadata().iloc[[0]].copy()
        metadata.loc[0, "Accession"] = "001"

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = self.write_fasta(folder, ">001\nACGT\n")
            metadata_path = self.write_metadata(folder, metadata, ".csv")
            data, errors, _ = load_and_validate_inputs(fasta_path, metadata_path)

        self.assertEqual(errors, [])
        self.assertEqual(data.loc[0, "Accession"], "001")
        self.assertEqual(data.loc[0, "Sample_ID"], "001")

    def test_ncbi_metadata_aliases_are_standardized(self) -> None:
        metadata = self.valid_metadata().rename(
            columns={
                "Pathogen": "Organism_Name",
                "Gene_Segment": "Segment",
                "Location_Name": "Geo_Location",
                "Host_Source": "Host",
            }
        )

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = self.write_fasta(
                folder, ">ACC001\nACGT\n>ACC002\nACGT\n"
            )
            metadata_path = self.write_metadata(folder, metadata, ".csv")
            data, errors, _ = load_and_validate_inputs(fasta_path, metadata_path)

        self.assertEqual(errors, [])
        self.assertEqual(data["Pathogen"].tolist(), ["Pathogen A", "Pathogen B"])
        self.assertEqual(data["Gene_Segment"].tolist(), ["Segment 1", "Complete genome"])
        self.assertEqual(data["Location_Name"].tolist(), ["Location A", "Location B"])
        self.assertEqual(data["Host_Source"].tolist(), ["Human", "Environmental"])

    def test_accession_version_matches_unversioned_metadata(self) -> None:
        metadata = self.valid_metadata().iloc[[0]].copy()
        metadata.loc[0, "Accession"] = "QB061804"

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = self.write_fasta(folder, ">QB061804.1\nACGT\n")
            metadata_path = self.write_metadata(folder, metadata, ".csv")
            data, errors, _ = load_and_validate_inputs(fasta_path, metadata_path)

        self.assertEqual(errors, [])
        self.assertEqual(data.loc[0, "Accession"], "QB061804")
        self.assertEqual(data.loc[0, "Sample_ID"], "QB061804.1")

    def test_missing_coordinates_are_allowed_with_geo_location(self) -> None:
        metadata = self.valid_metadata().drop(columns=["Latitude", "Longitude"])
        metadata = metadata.rename(
            columns={
                "Pathogen": "Organism_Name",
                "Gene_Segment": "Segment",
                "Location_Name": "Geo_Location",
                "Host_Source": "Host",
            }
        )

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = self.write_fasta(
                folder, ">ACC001.1\nACGT\n>ACC002.1\nACGT\n"
            )
            metadata_path = self.write_metadata(folder, metadata, ".csv")
            data, errors, _ = load_and_validate_inputs(fasta_path, metadata_path)

        self.assertEqual(errors, [])
        self.assertTrue(data["Latitude"].isna().all())
        self.assertTrue(data["Longitude"].isna().all())
        self.assertEqual(data["Location_Name"].tolist(), ["Location A", "Location B"])

    def test_unmatched_accessions_are_reported_in_both_directions(self) -> None:
        metadata = self.valid_metadata()
        metadata.loc[1, "Accession"] = "META_ONLY"

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = self.write_fasta(
                folder, ">ACC001\nACGT\n>FASTA_ONLY\nACGT\n"
            )
            metadata_path = self.write_metadata(folder, metadata, ".csv")
            data, errors, _ = load_and_validate_inputs(fasta_path, metadata_path)

        combined_errors = "\n".join(errors)
        self.assertEqual(len(data), 2)
        self.assertIn("META_ONLY", combined_errors)
        self.assertIn("FASTA_ONLY", combined_errors)
        self.assertIn("without a matching FASTA", combined_errors)
        self.assertIn("without matching metadata", combined_errors)

    def test_duplicate_accessions_are_reported(self) -> None:
        metadata = self.valid_metadata()
        metadata.loc[1, "Accession"] = "ACC001.2"

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = self.write_fasta(
                folder, ">ACC001.1\nACGT\n>ACC001.2\nACGT\n"
            )
            metadata_path = self.write_metadata(folder, metadata, ".xlsx")
            _, errors, _ = load_and_validate_inputs(fasta_path, metadata_path)

        combined_errors = "\n".join(errors)
        self.assertIn("Duplicate metadata accession", combined_errors)
        self.assertIn("Duplicate FASTA accession", combined_errors)

    def test_invalid_fasta_sequence_is_reported(self) -> None:
        metadata = self.valid_metadata().iloc[[0]].copy()

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = self.write_fasta(folder, ">ACC001\nACGX\n")
            metadata_path = self.write_metadata(folder, metadata, ".csv")
            _, errors, _ = load_and_validate_inputs(fasta_path, metadata_path)

        self.assertIn("Invalid sequence character(s) 'X'", "\n".join(errors))

    def test_bad_metadata_values_are_reported_without_deleting_rows(self) -> None:
        metadata = self.valid_metadata()
        metadata.loc[0, "Pathogen"] = " "
        metadata.loc[1, "Collection_Date"] = "not a date"
        metadata.loc[0, "Latitude"] = 95
        metadata.loc[1, "Longitude"] = "west"

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = self.write_fasta(
                folder, ">ACC001\nACGT\n>ACC002\nACGT\n"
            )
            metadata_path = self.write_metadata(folder, metadata, ".xlsx")
            data, errors, _ = load_and_validate_inputs(fasta_path, metadata_path)

        combined_errors = "\n".join(errors)
        self.assertEqual(len(data), 2)
        self.assertIn("Missing Pathogen", combined_errors)
        self.assertIn("Invalid Collection_Date", combined_errors)
        self.assertIn("Latitude outside", combined_errors)
        self.assertIn("Non-numeric Longitude", combined_errors)

    def test_ncbi_numeric_segment_with_spaces_and_mixed_case(self) -> None:
        metadata = self.valid_metadata().iloc[[0]].copy()
        metadata = metadata.rename(
            columns={
                "Pathogen": "  ORGANISM NAME  ",
                "Gene_Segment": "  sEgMeNt  ",
                "Location_Name": " Geo Location ",
                "Host_Source": " hOsT ",
            }
        )
        metadata.iloc[0, metadata.columns.get_loc("  sEgMeNt  ")] = "4"

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = self.write_fasta(folder, ">ACC001 segment 4 HA gene\nACGT\n")
            metadata_path = self.write_metadata(folder, metadata, ".csv")
            data, errors, _ = load_and_validate_inputs(fasta_path, metadata_path)

        self.assertEqual(errors, [])
        self.assertEqual(data.loc[0, "Gene_Segment"], "Segment 4")
        self.assertEqual(data.loc[0, "Gene_Segment_Source"], "Metadata")

    def test_gene_and_region_aliases_are_accepted(self) -> None:
        examples = [
            (" Gene Name ", "HA", "Gene HA"),
            ("gEnE", "NA", "Gene NA"),
            (" Region ", "complete genome", "Complete genome"),
            (" gene_segment ", "Segment 4", "Segment 4"),
        ]
        for heading, value, expected in examples:
            with self.subTest(heading=heading):
                metadata = self.valid_metadata().iloc[[0]].copy()
                metadata = metadata.rename(columns={"Gene_Segment": heading})
                metadata.iloc[0, metadata.columns.get_loc(heading)] = value
                with tempfile.TemporaryDirectory() as folder:
                    fasta_path = self.write_fasta(folder, ">ACC001\nACGT\n")
                    metadata_path = self.write_metadata(folder, metadata, ".csv")
                    data, errors, _ = load_and_validate_inputs(
                        fasta_path, metadata_path
                    )
                self.assertEqual(errors, [])
                self.assertEqual(data.loc[0, "Gene_Segment"], expected)

    def test_fasta_header_supplies_segment_when_metadata_has_none(self) -> None:
        metadata = self.valid_metadata().iloc[[0]].drop(columns="Gene_Segment")

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = self.write_fasta(
                folder, ">ACC001 segment 4 hemagglutinin (HA) gene\nACGT\n"
            )
            metadata_path = self.write_metadata(folder, metadata, ".csv")
            data, errors, _ = load_and_validate_inputs(fasta_path, metadata_path)

        self.assertEqual(errors, [])
        self.assertEqual(data.loc[0, "Gene_Segment"], "Segment 4")
        self.assertEqual(data.loc[0, "Gene_Segment_Source"], "FASTA header")

    def test_header_gene_name_and_acronym_without_gene_word(self) -> None:
        metadata = self.valid_metadata().iloc[[0]].drop(columns="Gene_Segment")

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = self.write_fasta(folder, ">ACC001 hemagglutinin (HA)\nACGT\n")
            metadata_path = self.write_metadata(folder, metadata, ".csv")
            data, errors, _ = load_and_validate_inputs(fasta_path, metadata_path)

        self.assertEqual(errors, [])
        self.assertEqual(data.loc[0, "Gene_Segment"], "Gene HA")

    def test_complete_genome_header_for_nonsegmented_pathogen(self) -> None:
        metadata = self.valid_metadata().drop(columns="Gene_Segment")
        metadata["Pathogen"] = "One nonsegmented virus"

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = self.write_fasta(
                folder,
                ">ACC001 complete genome\nACGT\n>ACC002 complete genome\nACGA\n",
            )
            metadata_path = self.write_metadata(folder, metadata, ".csv")
            data, errors, _ = load_and_validate_inputs(fasta_path, metadata_path)

        self.assertEqual(errors, [])
        self.assertEqual(data["Gene_Segment"].tolist(), ["Complete genome"] * 2)
        self.assertEqual(data["Gene_Segment_Source"].tolist(), ["FASTA header"] * 2)

    def test_explicit_complete_genome_confirmation_labels_unknown_regions(self) -> None:
        metadata = self.valid_metadata().drop(columns="Gene_Segment")
        metadata["Pathogen"] = "One pathogen"

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = self.write_fasta(folder, ">ACC001\nACGT\n>ACC002\nACGA\n")
            metadata_path = self.write_metadata(folder, metadata, ".csv")
            data, errors, _ = load_and_validate_inputs(
                fasta_path, metadata_path, complete_genomes=True
            )

        self.assertEqual(errors, [])
        self.assertEqual(data["Gene_Segment"].tolist(), ["Complete genome"] * 2)
        self.assertEqual(
            data["Gene_Segment_Source"].tolist(),
            ["User-confirmed complete genome"] * 2,
        )

    def test_complete_genome_mode_rejects_segment_annotations(self) -> None:
        metadata = self.valid_metadata().iloc[[0]].copy()
        metadata["Pathogen"] = "One pathogen"

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = self.write_fasta(folder, ">ACC001 segment 1\nACGT\n")
            metadata_path = self.write_metadata(folder, metadata, ".csv")
            _, errors, _ = load_and_validate_inputs(
                fasta_path, metadata_path, complete_genomes=True
            )

        self.assertIn("Complete-genome mode conflicts", "\n".join(errors))

    def test_complete_genome_mode_rejects_fasta_only_segment(self) -> None:
        metadata = self.valid_metadata().iloc[[0]].drop(columns="Gene_Segment")

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = self.write_fasta(folder, ">ACC001 segment 1\nACGT\n")
            metadata_path = self.write_metadata(folder, metadata, ".csv")
            _, errors, _ = load_and_validate_inputs(
                fasta_path, metadata_path, complete_genomes=True
            )

        self.assertIn("conflicts with the FASTA header", "\n".join(errors))

    def test_complete_genome_mode_rejects_partial_sequence_evidence(self) -> None:
        metadata = self.valid_metadata().iloc[[0]].drop(columns="Gene_Segment")
        metadata["Nuc_Completeness"] = "partial"

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = self.write_fasta(
                folder, ">ACC001 complete cds\nACGT\n"
            )
            metadata_path = self.write_metadata(folder, metadata, ".csv")
            _, errors, _ = load_and_validate_inputs(
                fasta_path, metadata_path, complete_genomes=True
            )

        error_text = "\n".join(errors)
        self.assertIn("Nuc_Completeness", error_text)
        self.assertIn("partial sequence or CDS", error_text)

    def test_complete_genome_mode_requires_one_pathogen(self) -> None:
        metadata = self.valid_metadata().drop(columns="Gene_Segment")

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = self.write_fasta(folder, ">ACC001\nACGT\n>ACC002\nACGA\n")
            metadata_path = self.write_metadata(folder, metadata, ".csv")
            _, errors, _ = load_and_validate_inputs(
                fasta_path, metadata_path, complete_genomes=True
            )

        self.assertIn("requires all samples to have one pathogen", "\n".join(errors))

    def test_unannotated_single_pathogen_requires_confirmation(self) -> None:
        metadata = self.valid_metadata().copy().drop(columns="Gene_Segment")
        metadata["Pathogen"] = "One pathogen"

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = self.write_fasta(folder, ">ACC001\nACGT\n>ACC002\nACGTA\n")
            metadata_path = self.write_metadata(folder, metadata, ".csv")
            _, errors_without, _ = load_and_validate_inputs(fasta_path, metadata_path)
            data, errors_with, _ = load_and_validate_inputs(
                fasta_path, metadata_path, confirm_same_region=True
            )

        self.assertIn("Cannot determine genomic region", "\n".join(errors_without))
        self.assertEqual(errors_with, [])
        self.assertEqual(data["Gene_Segment"].tolist(), ["Unspecified_Region"] * 2)
        self.assertEqual(
            data["Gene_Segment_Source"].tolist(),
            ["User-confirmed unspecified region"] * 2,
        )

    def test_default_region_used_only_when_no_explicit_annotation(self) -> None:
        metadata = self.valid_metadata().copy().drop(columns="Gene_Segment")
        metadata["Pathogen"] = "One pathogen"

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = self.write_fasta(
                folder, ">ACC001 segment 4\nACGT\n>ACC002\nACGA\n"
            )
            metadata_path = self.write_metadata(folder, metadata, ".csv")
            data, errors, _ = load_and_validate_inputs(
                fasta_path, metadata_path, default_gene_segment="Segment 4"
            )

        self.assertEqual(errors, [])
        self.assertEqual(data["Gene_Segment"].tolist(), ["Segment 4"] * 2)
        self.assertEqual(
            data["Gene_Segment_Source"].tolist(), ["FASTA header", "User override"]
        )

    def test_conflicting_metadata_and_header_are_rejected(self) -> None:
        metadata = self.valid_metadata().iloc[[0]].copy()
        metadata.loc[0, "Gene_Segment"] = "Segment 4"

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = self.write_fasta(folder, ">ACC001 segment 5\nACGT\n")
            metadata_path = self.write_metadata(folder, metadata, ".csv")
            _, errors, _ = load_and_validate_inputs(fasta_path, metadata_path)

        self.assertIn("Metadata and FASTA header disagree", "\n".join(errors))

    def test_secondary_gene_metadata_cannot_contradict_fasta(self) -> None:
        metadata = self.valid_metadata().iloc[[0]].copy()
        metadata.loc[0, "Gene_Segment"] = "Segment 4"
        metadata["Gene_Name"] = "HA"

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = self.write_fasta(
                folder, ">ACC001 segment 4 neuraminidase (NA) gene\nACGT\n"
            )
            metadata_path = self.write_metadata(folder, metadata, ".csv")
            _, errors, _ = load_and_validate_inputs(fasta_path, metadata_path)

        self.assertIn("Metadata and FASTA header disagree", "\n".join(errors))

    def test_confirmation_does_not_merge_multiple_pathogens(self) -> None:
        metadata = self.valid_metadata().drop(columns="Gene_Segment")

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = self.write_fasta(folder, ">ACC001\nACGT\n>ACC002\nACGT\n")
            metadata_path = self.write_metadata(folder, metadata, ".csv")
            _, errors, _ = load_and_validate_inputs(
                fasta_path, metadata_path, confirm_same_region=True
            )

        self.assertIn("one pathogen only", "\n".join(errors))

    def test_unspecified_region_cannot_be_used_without_confirmation(self) -> None:
        metadata = self.valid_metadata().iloc[[0]].copy()
        metadata.loc[0, "Gene_Segment"] = "Unspecified_Region"

        with tempfile.TemporaryDirectory() as folder:
            fasta_path = self.write_fasta(folder, ">ACC001\nACGT\n")
            metadata_path = self.write_metadata(folder, metadata, ".csv")
            _, errors, _ = load_and_validate_inputs(fasta_path, metadata_path)

        self.assertIn("requires --confirm-same-region", "\n".join(errors))


if __name__ == "__main__":
    unittest.main()
