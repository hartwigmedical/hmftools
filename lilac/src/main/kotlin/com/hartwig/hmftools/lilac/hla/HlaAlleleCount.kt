package com.hartwig.hmftools.lilac.hla

import com.hartwig.hmftools.lilac.read.FragmentSequences

data class HlaAlleleCount(val allele: HlaAllele, val uniqueCoverage: Int, val combinedCoverage: Double): Comparable<HlaAlleleCount> {

    companion object {

        fun proteinCoverage(fragmentSequences: List<FragmentSequences>): List<HlaAlleleCount> {
            return create(fragmentSequences) {it}
        }

        fun groupCoverage(fragmentSequences: List<FragmentSequences>): List<HlaAlleleCount> {
            return create(fragmentSequences) { it.alleleGroup()}
        }

        fun create(fragmentSequences: List<FragmentSequences>, type: (HlaAllele) -> HlaAllele): List<HlaAlleleCount> {
            val result = mutableListOf<HlaAlleleCount>()

            val uniqueCoverageMap = mutableMapOf<HlaAllele, Int>()
            val combinedCoverageMap = mutableMapOf<HlaAllele, Double>()

            // Counts
            for (fragment in fragmentSequences) {
                val fullAlleles = fragment.full.map(type).toSet()
                val partialAlleles = fragment.partial.map(type).toSet()

                if (fullAlleles.size == 1 && partialAlleles.isEmpty())  {
                    uniqueCoverageMap.compute(fullAlleles.first()) {_, oldValue ->  (oldValue ?: 0) + 1}
                } else {
                    val contribution = 1.0 / (fullAlleles.size + partialAlleles.size)
                    fullAlleles.forEach {combinedCoverageMap.compute(it) {_, oldValue -> (oldValue ?: 0.0) + contribution} }
                    partialAlleles.forEach {combinedCoverageMap.compute(it) {_, oldValue -> (oldValue ?: 0.0) + contribution} }
                }
            }

            // Produce results
            val hlaAlleles = uniqueCoverageMap.keys union combinedCoverageMap.keys
            for (allele in hlaAlleles) {
                val uniqueCoverage = uniqueCoverageMap[allele] ?: 0
                val combinedCoverage = combinedCoverageMap[allele] ?: 0.0

                result.add(HlaAlleleCount(allele, uniqueCoverage, combinedCoverage))
            }

            return result.sortedDescending()
        }
    }

    override fun compareTo(other: HlaAlleleCount): Int {
        val uniqueCompare = uniqueCoverage.compareTo(other.uniqueCoverage)
        if (uniqueCompare != 0) {
            return uniqueCompare
        }
        return combinedCoverage.compareTo(other.combinedCoverage)
    }
}