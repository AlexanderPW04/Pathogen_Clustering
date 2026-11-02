"""Population structure of pooled samples: Yue-Clayton dissimilarity and PCoA.

Paper method (Fontenele et al. 2021, Section 2.6), reproduced here:

1. Represent each sample by A/C/G/T frequencies at every reference position.
2. At each position, compute the Yue-Clayton dissimilarity between two
   samples' four-base frequency vectors.
3. If a frequency cannot be calculated at a position for either sample (zero
   depth or the position is absent from the table), that position contributes
   0 - the paper's stated treatment of zero-depth sites.
4. Sum the per-position values to get the pairwise distance, then run a
   principal coordinates analysis (PCoA; Gower 1966).

The paper does not print the Yue-Clayton formula. This module uses the
standard definition (Yue & Clayton 2005):

    similarity    = sum(a*b) / (sum(a^2) + sum(b^2) - sum(a*b))
    dissimilarity = 1 - similarity

This needs per-position base counts. A table of called SNVs alone does not
contain reference-allele frequencies at unreported positions, so the exact
PCoA is disabled for SNV-only input and nothing is invented to replace it.
"""

import numpy as np
import pandas as pd

from genomics.pooled_observations import COUNT_COLUMNS


PAPER_PCOA = "Paper method: Yue-Clayton PCoA of A/C/G/T frequencies"
EXPLORATORY_ORDINATION = "Exploratory: PCoA of Jaccard distances between reported SNV sets"


def yue_clayton_dissimilarity(frequencies_a, frequencies_b):
    """Yue-Clayton dissimilarity of two frequency vectors (0 = identical)."""

    a = np.asarray(frequencies_a, dtype=float)
    b = np.asarray(frequencies_b, dtype=float)
    if np.sum(a) <= 0 or np.sum(b) <= 0:
        return None  # a frequency cannot be calculated without reads
    shared = float(np.sum(a * b))
    denominator = float(np.sum(a * a) + np.sum(b * b) - shared)
    return 1.0 - shared / denominator


def frequency_arrays(position_table, sample_ids, reference_id):
    """Build an array [sample, position, base] of frequencies; NaN = no reads."""

    rows = position_table[
        position_table["Reference_ID"].eq(reference_id)
        & position_table["Sample_ID"].isin(sample_ids)
    ]
    last_position = int(rows["Position"].max()) if not rows.empty else 0
    frequencies = np.full((len(sample_ids), last_position, 4), np.nan, dtype=np.float32)
    index_by_sample = {sample_id: index for index, sample_id in enumerate(sample_ids)}

    for sample_id, group in rows.groupby("Sample_ID"):
        counts = group[COUNT_COLUMNS].to_numpy(dtype=np.float64)
        totals = counts.sum(axis=1)
        positions = group["Position"].to_numpy() - 1
        has_reads = totals > 0
        values = np.full(counts.shape, np.nan)
        values[has_reads] = counts[has_reads] / totals[has_reads, None]
        frequencies[index_by_sample[sample_id], positions] = values
    return frequencies


def paper_distance_matrix(frequencies):
    """Sum per-position Yue-Clayton dissimilarities; uncomputable positions = 0."""

    sample_count = frequencies.shape[0]
    usable = ~np.isnan(frequencies[:, :, 0])
    squares = np.where(usable, np.nansum(frequencies ** 2, axis=2), 0.0)
    distances = np.zeros((sample_count, sample_count))
    compared_positions = np.zeros((sample_count, sample_count), dtype=int)

    for first in range(sample_count):
        for second in range(first + 1, sample_count):
            both = usable[first] & usable[second]
            shared = np.sum(frequencies[first, both] * frequencies[second, both], axis=1)
            denominator = squares[first, both] + squares[second, both] - shared
            per_position = 1.0 - shared / denominator
            total = float(np.sum(per_position))
            distances[first, second] = distances[second, first] = total
            compared_positions[first, second] = compared_positions[second, first] = int(both.sum())
    return distances, compared_positions


def jaccard_distance_matrix(variant_sets):
    """1 - |shared| / |either| between reported variant sets."""

    sample_count = len(variant_sets)
    distances = np.zeros((sample_count, sample_count))
    for first in range(sample_count):
        for second in range(first + 1, sample_count):
            either = variant_sets[first] | variant_sets[second]
            shared = variant_sets[first] & variant_sets[second]
            value = 1.0 - len(shared) / len(either) if either else 0.0
            distances[first, second] = distances[second, first] = value
    return distances


def principal_coordinates(distances):
    """Classical PCoA (Gower 1966) with an explicit negative-eigenvalue report.

    B = -1/2 * J D^2 J is decomposed; axes come from positive eigenvalues.
    Non-Euclidean distances can give negative eigenvalues; they are reported,
    not corrected.
    """

    distances = np.asarray(distances, dtype=float)
    count = distances.shape[0]
    centering = np.eye(count) - np.ones((count, count)) / count
    gower = -0.5 * centering @ (distances ** 2) @ centering
    eigenvalues, eigenvectors = np.linalg.eigh(gower)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]

    tolerance = 1e-9 * max(1.0, float(np.max(np.abs(eigenvalues))))
    positive = eigenvalues > tolerance
    negative = eigenvalues < -tolerance
    coordinates = eigenvectors[:, positive] * np.sqrt(eigenvalues[positive])
    positive_sum = float(np.sum(eigenvalues[positive]))
    explained = eigenvalues[positive] / positive_sum if positive_sum > 0 else eigenvalues[positive]
    return {
        "coordinates": coordinates,
        "eigenvalues": eigenvalues,
        "explained": explained,
        "negative_eigenvalues": int(negative.sum()),
        "negative_magnitude_ratio": (
            float(np.sum(np.abs(eigenvalues[negative])) / positive_sum) if positive_sum > 0 else 0.0
        ),
    }


def nucleotide_diversity(frequencies):
    """Mean per-position diversity 1 - sum(f^2) over positions with reads."""

    usable = ~np.isnan(frequencies[:, :, 0])
    values = []
    for index in range(frequencies.shape[0]):
        rows = frequencies[index, usable[index]]
        values.append(float(np.mean(1.0 - np.sum(rows ** 2, axis=1))) if len(rows) else np.nan)
    return values


def coordinates_table(sample_info, pcoa):
    """Combine sample metadata with the first two PCoA axes."""

    table = sample_info[["Sample_ID", "Collection_Date", "Location_Name"]].copy().reset_index(drop=True)
    coordinates = pcoa["coordinates"]
    for axis in range(2):
        table[f"PCoA_{axis + 1}"] = coordinates[:, axis] if coordinates.shape[1] > axis else 0.0
    return table


def population_structure(samples, variants, position_table=None, minimum_samples=3):
    """Run the paper's PCoA when possible, else a labeled exploratory ordination."""

    included = samples[samples["Analysis_Included"]].sort_values(["Collection_Date", "Sample_ID"])
    result = {"method": None, "available": False, "message": "", "pcoa": None,
              "coordinates": None, "distances": None, "diversity": None}

    # Positions from different references cannot share one matrix, so only the
    # reference used by the most samples is analyzed; the others are named.
    single_reference = included[~included["Reference_ID"].astype(str).str.contains("; ")]
    reference_counts = single_reference["Reference_ID"].value_counts()
    if reference_counts.empty:
        result["message"] = "No pooled sample has variant records on a single reference."
        return result
    reference_id = reference_counts.index[0]
    left_out = sorted(set(included["Sample_ID"]) - set(
        single_reference.loc[single_reference["Reference_ID"].eq(reference_id), "Sample_ID"]
    ))
    included = included[included["Reference_ID"].eq(reference_id)]
    result["reference_id"] = reference_id
    result["left_out"] = left_out

    if len(included) < minimum_samples:
        result["message"] = (
            f"Population structure needs at least {minimum_samples} pooled samples that "
            f"passed quality control on one reference; {len(included)} are available."
        )
        return result

    sample_ids = list(included["Sample_ID"])
    has_counts = position_table is not None and position_table.attrs.get("has_base_counts", False)
    counted = set()
    if has_counts:
        counted = set(position_table.loc[position_table["Reference_ID"].eq(reference_id), "Sample_ID"])

    if has_counts and set(sample_ids) <= counted:
        frequencies = frequency_arrays(position_table, sample_ids, reference_id)
        distances, compared = paper_distance_matrix(frequencies)
        result["method"] = PAPER_PCOA
        result["diversity"] = pd.DataFrame({
            "Sample_ID": sample_ids,
            "Mean_Nucleotide_Diversity": nucleotide_diversity(frequencies),
            "Positions_With_Reads": (~np.isnan(frequencies[:, :, 0])).sum(axis=1),
        })
        result["message"] = (
            "Computed with the paper's method from supplied A/C/G/T counts. Positions "
            "without reads in either sample contribute 0, as stated in the paper."
        )
        result["compared_positions"] = compared
    else:
        included_variants = variants[variants["Analysis_Included"] & variants["Sample_ID"].isin(sample_ids)]
        variant_sets = [
            set(included_variants.loc[included_variants["Sample_ID"].eq(sample_id), "Variant_Key"])
            for sample_id in sample_ids
        ]
        distances = jaccard_distance_matrix(variant_sets)
        result["method"] = EXPLORATORY_ORDINATION
        missing = "no per-position A/C/G/T counts were supplied"
        if has_counts:
            missing = "A/C/G/T counts were not supplied for every included sample"
        result["message"] = (
            f"The paper's Yue-Clayton PCoA cannot be computed because {missing}; a table of "
            "called SNVs does not contain reference-allele frequencies elsewhere. This "
            "exploratory ordination compares reported SNV sets only. An unreported SNV "
            "is not evidence of absence, so sequencing depth can drive these distances."
        )

    if left_out:
        result["message"] += (
            f" Left out because they use a different reference than {reference_id}: "
            + ", ".join(left_out) + "."
        )
    pcoa = principal_coordinates(distances)
    result["available"] = True
    result["pcoa"] = pcoa
    result["coordinates"] = coordinates_table(included, pcoa)
    result["distances"] = pd.DataFrame(distances, index=sample_ids, columns=sample_ids)
    return result
