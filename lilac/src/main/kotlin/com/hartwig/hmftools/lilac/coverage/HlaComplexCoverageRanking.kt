package com.hartwig.hmftools.lilac.coverage

import com.hartwig.hmftools.lilac.hla.HlaAllele
import kotlin.math.min

/**
 * Complexes within [maxDistanceFromTopScore] fragments of the highest aligned fragment count are considered as possible solutions.
 * These are then ranked by:
 * - first prioritising solutions with the fewest alleles with wildcard matches,
 * - then choosing solutions with the most homozygous alleles,
 * - then choosing solution with the most common alleles,
 * - then choosing solution with the least recovered alleles
 * - finally choosing the solution with the lowest number.
 */
class HlaComplexCoverageRanking(private val maxDistanceFromTopScore: Int, private val common: List<HlaAllele>, private val recovered: List<HlaAllele>) {

    fun candidateRanking(complexes: List<HlaComplexCoverage>): List<HlaComplexCoverage> {
        require(complexes.isNotEmpty())

        val topScore = complexes.map { it.totalCoverage }.max()!!
        if (topScore == 0) {
            return listOf()
        }

        return complexes
                .filter { it.totalCoverage >= topScore - maxDistanceFromTopScore }
                .sortedWith(Comparator { x, y -> compare(x, y) })
    }

    private fun compare(o1: HlaComplexCoverage, o2: HlaComplexCoverage): Int {
        val wildcardCount = o1.wildcardCount().compareTo(o2.wildcardCount())
        if (wildcardCount != 0) {
            return wildcardCount
        }

        val homozygousCompare = o1.homozygousAlleles().compareTo(o2.homozygousAlleles())
        if (homozygousCompare != 0) {
            return -homozygousCompare
        }

        val o1CommonCount = o1.commonCount()
        val o2CommonCount = o2.commonCount()
        val commonCountCompare = o1CommonCount.compareTo(o2CommonCount)
        if (commonCountCompare != 0) {
            return -commonCountCompare
        }

        val o1RecoveredCount = o1.recoveredCount()
        val o2RecoveredCount = o2.recoveredCount()
        val recoveredCountCompare = o1RecoveredCount.compareTo(o2RecoveredCount)
        if (recoveredCountCompare != 0) {
            return recoveredCountCompare
        }

        for (i in 0 until min(o1.alleleCoverage.size, o2.alleleCoverage.size)) {
            val o1Allele = o1.alleleCoverage[i].allele
            val o2Allele = o2.alleleCoverage[i].allele
            val alleleCompare = o1Allele.compareTo(o2Allele)
            if (alleleCompare != 0) {
                return alleleCompare
            }
        }

        throw UnsupportedOperationException("Should not be able to make it to here")
    }

    private fun HlaComplexCoverage.wildcardCount(): Int {
        return this.alleleCoverage.filter { it.wildCoverage > 0 }.count()
    }

    private fun HlaComplexCoverage.commonCount(): Int {
        return this.alleleCoverage.map { it.allele }.filter { it in common }.count()
    }

    private fun HlaComplexCoverage.recoveredCount(): Int {
        return this.alleleCoverage.map { it.allele }.filter { it in recovered }.count()
    }

}