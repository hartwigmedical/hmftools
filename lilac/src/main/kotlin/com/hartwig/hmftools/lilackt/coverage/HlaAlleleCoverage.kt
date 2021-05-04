package com.hartwig.hmftools.lilackt.coverage

import com.hartwig.hmftools.lilackt.hla.HlaAllele
import com.hartwig.hmftools.lilackt.read.FragmentAlleles
import kotlin.math.roundToInt

data class HlaAlleleCoverage(val allele: HlaAllele, val uniqueCoverage: Int, val sharedCoverage: Double, val wildCoverage: Double) : Comparable<HlaAlleleCoverage> {
    val totalCoverage = uniqueCoverage + sharedCoverage + wildCoverage

    companion object {

        fun List<HlaAlleleCoverage>.expand(): List<HlaAlleleCoverage> {
            val a = this.filter { it.allele.gene == "A" }
            val b = this.filter { it.allele.gene == "B" }
            val c = this.filter { it.allele.gene == "C" }

            return (a.splitSingle() + b.splitSingle() + c.splitSingle()).sortedBy { it.allele }
        }

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

                if (fullAlleles.size == 1 && partialAlleles.isEmpty())  {
                    uniqueCoverageMap.compute(fullAlleles.first()) {_, oldValue ->  (oldValue ?: 0) + 1}
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

        private fun List<HlaAlleleCoverage>.splitSingle(): List<HlaAlleleCoverage> {
            if (this.size == 1) {
                val single = this[0]
                val first = HlaAlleleCoverage(single.allele, single.uniqueCoverage / 2, single.sharedCoverage / 2, single.wildCoverage / 2)
                val remainder = HlaAlleleCoverage(single.allele, single.uniqueCoverage - first.uniqueCoverage, single.sharedCoverage - first.sharedCoverage, single.wildCoverage - single.wildCoverage)
                return listOf(first, remainder).sortedBy { it.totalCoverage }.reversed()

            }

            return this
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