package com.hartwig.hmftools.lilac.coverage

import com.hartwig.hmftools.lilac.hla.HlaAllele
import junit.framework.Assert.assertEquals
import org.junit.Test


class HlaComplexCoverageRankingTest {

    val victim = HlaComplexCoverageRanking(3)

    val a1 = HlaAlleleCoverage(HlaAllele("A*01:01"), 10, 0.0, 0.0)
    val a2 = HlaAlleleCoverage(HlaAllele("A*01:02"), 10, 0.0, 0.0)
    val a3 = HlaAlleleCoverage(HlaAllele("A*01:03"), 10, 0.0, 0.0)
    val b1 = HlaAlleleCoverage(HlaAllele("B*01:01"), 10, 0.0, 0.0)
    val b2 = HlaAlleleCoverage(HlaAllele("B*01:02"), 10, 0.0, 0.0)
    val b3 = HlaAlleleCoverage(HlaAllele("B*01:03"), 10, 0.0, 0.0)
    val c1 = HlaAlleleCoverage(HlaAllele("C*01:01"), 10, 0.0, 0.0)
    val c2 = HlaAlleleCoverage(HlaAllele("C*01:02"), 10, 0.0, 0.0)
    val c3 = HlaAlleleCoverage(HlaAllele("C*01:03"), 10, 0.0, 0.0)

    @Test
    fun testHighestAllele() {
        val highest = HlaComplexCoverage.create(listOf(a1, a2, b1, b2, c1, c2))
        val lowerList = mutableListOf<HlaComplexCoverage>()

        lowerList.add(HlaComplexCoverage.create(listOf(a1, a3, b1, b2, c1, c2)))
        lowerList.add(HlaComplexCoverage.create(listOf(a1, a2, b1, b3, c1, c2)))
        lowerList.add(HlaComplexCoverage.create(listOf(a1, a2, b1, b2, c1, c3)))

        for (lower in lowerList) {
            val complexes = mutableListOf(highest, lower).shuffled()
            val winner = victim.candidateRanking(complexes)[0]
            assertEquals(highest, winner)
        }
    }

    @Test
    fun testHomozygousWinsIfAllElseEqual() {
        val het = HlaComplexCoverage.create(listOf(a1, a2, b1, b2, c1, c2))
        val homList = mutableListOf<HlaComplexCoverage>()

        homList.add(HlaComplexCoverage.create(listOf(a1, a1, b1, b2, c1, c2)))
        homList.add(HlaComplexCoverage.create(listOf(a2, a2, b1, b2, c1, c2)))
        homList.add(HlaComplexCoverage.create(listOf(a1, a2, b1, b1, c1, c2)))
        homList.add(HlaComplexCoverage.create(listOf(a1, a2, b2, b2, c1, c2)))
        homList.add(HlaComplexCoverage.create(listOf(a1, a2, b1, b2, c1, c1)))
        homList.add(HlaComplexCoverage.create(listOf(a1, a2, b1, b2, c2, c2)))

        for (hom in homList) {
            val complexes = mutableListOf(het, hom).shuffled()
            val winner = victim.candidateRanking(complexes)[0]
            assertEquals(hom, winner)
        }

        val homA = HlaComplexCoverage.create(listOf(a1, a1, b1, b2, c1, c2))
        val homAAndB = HlaComplexCoverage.create(listOf(a1, a1, b2, b2, c1, c2))
        val complexes = mutableListOf(homA, homAAndB).shuffled()
        assertEquals(homAAndB, victim.candidateRanking(complexes)[0])
    }

    @Test
    fun testWildCoverage() {
        val het = HlaComplexCoverage.create(listOf(a1, a2, b1, b2, c1, c2))
        val homWithoutWild = HlaComplexCoverage.create(listOf(a1, a1, b1, b1, c1, c1))
        val homWithWild = HlaComplexCoverage.create(listOf(a1, a1, b1, b1, c1, c1.copy(uniqueCoverage = 9, wildCoverage = 1.0)))

        assertEquals(homWithoutWild, victim.candidateRanking(mutableListOf(het, homWithoutWild).shuffled())[0])
        assertEquals(het, victim.candidateRanking(mutableListOf(het, homWithWild).shuffled())[0])
    }


}
