package com.hartwig.hmftools.lilac.hla

import com.hartwig.hmftools.lilac.read.FragmentAlleles
import kotlin.math.roundToInt

data class HlaAlleleCoverage(val allele: HlaAllele, val uniqueCoverage: Int, val sharedCoverage: Double) : Comparable<HlaAlleleCoverage> {

    companion object {

        fun proteinCoverage(fragmentSequences: List<FragmentAlleles>): List<HlaAlleleCoverage> {
            return create(fragmentSequences) { it }
        }

        fun groupCoverage(fragmentSequences: List<FragmentAlleles>): List<HlaAlleleCoverage> {
            return create(fragmentSequences) { it.alleleGroup() }
        }

        fun create(fragmentSequences: List<FragmentAlleles>, type: (HlaAllele) -> HlaAllele): List<HlaAlleleCoverage> {
            val result = mutableListOf<HlaAlleleCoverage>()

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

                result.add(HlaAlleleCoverage(allele, uniqueCoverage, combinedCoverage))
            }

            return result.sortedDescending()
        }
    }

    override fun compareTo(other: HlaAlleleCoverage): Int {
        val uniqueCompare = uniqueCoverage.compareTo(other.uniqueCoverage)
        if (uniqueCompare != 0) {
            return uniqueCompare
        }
        return sharedCoverage.compareTo(other.sharedCoverage)
    }

    override fun toString(): String {
        return "$allele:[t=${(uniqueCoverage + sharedCoverage).roundToInt()}, u=$uniqueCoverage, s=${sharedCoverage.roundToInt()}]"
    }


}