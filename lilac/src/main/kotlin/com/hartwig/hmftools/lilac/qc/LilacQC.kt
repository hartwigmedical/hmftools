package com.hartwig.hmftools.lilac.qc

import org.apache.logging.log4j.LogManager


data class LilacQC(val aTypes: Int, val bTypes: Int, val cTypes: Int,
                   val unmatchedFragments: Int, val uniqueFragments: Int, val sharedFragments: Int, val wildcardFragments: Int,
                   val discardedIndelFragments: Int, val discardedIndelMaxCount: Int,
                   val unmatchedAminoAcidSequences: Int, val unmatchedAminoAcidMaxCount: Int) {
    val fittedFragments = uniqueFragments + sharedFragments + wildcardFragments;
    val percentUnmatchedFragments = 1.0 * unmatchedFragments / fittedFragments

    companion object {
        val logger = LogManager.getLogger(this::class.java)
    }

}


// Was looking at the warnings again and was thinking we could write a QC file with the following rows:
//
//A type found (T/F)
//B type found (T/F)
//C type found (T/F)
//Total fragments
//% of unmatched fragments
//% fitted fragments
//% wildcard fragments
//% shared fragments
//% unique fragments
//Total discarded INDEL fragments
//Max discarded INDEL fragments
//Count unmatched Amino Acid
//Max unmatched Amino Acid count
//Count unmatched Haplotypes (>2)
//Max unmatched haplotype support
//# of dropped alleles with unique support
//Max support for dropped allele