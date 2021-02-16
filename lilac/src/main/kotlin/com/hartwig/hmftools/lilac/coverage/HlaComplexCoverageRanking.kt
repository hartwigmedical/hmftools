package com.hartwig.hmftools.lilac.coverage

import kotlin.math.min

/**
 * Complexes within [maxDistanceFromTopScore] fragments of the highest aligned fragment count are considered as possible solutions.
 * These are then ranked by:
 * - first prioritising solutions with the lowest number of wildcard matches,
 * - then choosing solutions with the most homozygous alleles,
 * - then choosing the solution with the highest allele population frequency.
 */
class HlaComplexCoverageRanking(private val maxDistanceFromTopScore: Int = 3) { // TODO: MOVE INTO CONFIG

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
        val wildCoverageCompare = o1.wildCoverage.compareTo(o2.wildCoverage)
        if (wildCoverageCompare != 0) {
            return wildCoverageCompare
        }

        val homozygousCompare = o1.homozygousAlleles().compareTo(o2.homozygousAlleles())
        if (homozygousCompare != 0) {
            return -homozygousCompare
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

}