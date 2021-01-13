package com.hartwig.hmftools.lilac.hla

import com.hartwig.hmftools.lilac.read.FragmentSequences

data class HlaAlleleCount(val allele: HlaAllele, val uniqueCoverage: Int, val combinedCoverage: Double) {

    companion object {

        fun proteinCoverage(fragmentSequences: List<FragmentSequences>): List<HlaAlleleCount> {
            return create(fragmentSequences) {it.alleles}
        }

        fun groupCoverage(fragmentSequences: List<FragmentSequences>): List<HlaAlleleCount> {
            return create(fragmentSequences) {it.alleleGroups}
        }

        fun create(fragmentSequences: List<FragmentSequences>, type: (FragmentSequences) -> Collection<HlaAllele>): List<HlaAlleleCount> {
            val result = mutableListOf<HlaAlleleCount>()

            val uniqueCoverageMap = mutableMapOf<HlaAllele, Int>()
            val combinedCoverageMap = mutableMapOf<HlaAllele, Double>()

            // Counts
            for (fragment in fragmentSequences) {
                val alleles = type(fragment)
                if (alleles.size == 1) {
                    uniqueCoverageMap.compute(alleles.first()) {_, oldValue ->  (oldValue ?: 0) + 1}
                } else {
                    val contribution = 1.0 / alleles.size
                    alleles.forEach {combinedCoverageMap.compute(it) {_, oldValue -> (oldValue ?: 0.0) + contribution} }
                }
            }

            // Produce results
            val hlaAlleles = uniqueCoverageMap.keys union combinedCoverageMap.keys
            for (allele in hlaAlleles) {
                val uniqueCoverage = uniqueCoverageMap[allele] ?: 0
                val combinedCoverage = combinedCoverageMap[allele] ?: 0.0

                result.add(HlaAlleleCount(allele, uniqueCoverage, combinedCoverage))
            }

            return result.sortedBy { -it.uniqueCoverage }
        }
    }
}