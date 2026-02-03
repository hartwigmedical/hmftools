package com.hartwig.hmftools.cider

// Determination of "somatic hypermutation status" from V-region sequence, which is a prognostic indicator for chronic lymphocytic leukemia.
// https://pmc.ncbi.nlm.nih.gov/articles/PMC7248390/

enum class SomaticHypermutationStatus {
    UNMUTATED,
    MUTATED_BORDERLINE,
    MUTATED;

    companion object
    {
        // Determine the status from the % identity of the V-region to the IMGT reference sequence.
        fun determineFromVGeneIdentity(pctIdentity: Double): SomaticHypermutationStatus
        {
            return if (pctIdentity >= 98) UNMUTATED
            else if (pctIdentity >= 97) MUTATED_BORDERLINE
            else MUTATED
        }
    }
}

// TODO: unit test!!!
