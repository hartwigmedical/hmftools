package com.hartwig.hmftools.lilac.hla

import com.hartwig.hmftools.lilac.read.FragmentAlleles
import kotlin.math.roundToInt

data class HlaAlleleCoverage(val allele: HlaAllele, val uniqueCoverage: Int, val sharedCoverage: Double, val wildCoverage: Double) : Comparable<HlaAlleleCoverage> {
    val totalCoverage = uniqueCoverage + sharedCoverage + wildCoverage

    companion object {

        fun proteinCoverage(fragmentSequences: List<FragmentAlleles>): List<HlaAlleleCoverage> {
            return create(fragmentSequences) { it }
        }

        fun groupCoverage(fragmentSequences: List<FragmentAlleles>): List<HlaAlleleCoverage> {
            return create(fragmentSequences) { it.asAlleleGroup() }
        }

        fun create(fragmentSequences: List<FragmentAlleles>, type: (HlaAllele) -> HlaAllele): List<HlaAlleleCoverage> {
            val result = mutableListOf<HlaAlleleCoverage>()

            val uniqueCoverageMap = mutableMapOf<HlaAllele, Int>()
            val combinedCoverageMap = mutableMapOf<HlaAllele, Double>()
            val wildCoverageMap = mutableMapOf<HlaAllele, Double>()

            // Counts
            for (fragment in fragmentSequences) {
                val fullAlleles = fragment.full.map(type).toSet()
                val partialAlleles = fragment.partial.map(type).toSet()
                val wildAlleles = fragment.wild.map(type).toSet()

                if (fullAlleles.size == 1 && partialAlleles.isEmpty() && wildAlleles.isEmpty())  {
                    uniqueCoverageMap.compute(fullAlleles.first()) {_, oldValue ->  (oldValue ?: 0) + 1}
//                    if (fullAlleles.first() == HlaAllele("C*07:57")) {
//                        println(fragment.aminoAcidFragment.id)
//                    }

                } else {
                    val contribution = 1.0 / (fullAlleles.size + partialAlleles.size + wildAlleles.size)
                    fullAlleles.forEach {combinedCoverageMap.compute(it) {_, oldValue -> (oldValue ?: 0.0) + contribution} }
                    partialAlleles.forEach {combinedCoverageMap.compute(it) {_, oldValue -> (oldValue ?: 0.0) + contribution} }
                    wildAlleles.forEach {wildCoverageMap.compute(it) {_, oldValue -> (oldValue ?: 0.0) + contribution} }
                }
            }

            // Produce results
            val hlaAlleles = uniqueCoverageMap.keys union combinedCoverageMap.keys
            for (allele in hlaAlleles) {
                val uniqueCoverage = uniqueCoverageMap[allele] ?: 0
                val combinedCoverage = combinedCoverageMap[allele] ?: 0.0
                val wildCoverage = wildCoverageMap[allele] ?: 0.0

                result.add(HlaAlleleCoverage(allele, uniqueCoverage, combinedCoverage, wildCoverage))
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
        return "$allele[${totalCoverage.roundToInt()},$uniqueCoverage,${sharedCoverage.roundToInt()},${wildCoverage.roundToInt()}]"
    }


}