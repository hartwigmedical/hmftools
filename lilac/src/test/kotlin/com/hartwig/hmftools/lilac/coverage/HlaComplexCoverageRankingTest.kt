package com.hartwig.hmftools.lilac.coverage

import com.hartwig.hmftools.lilac.hla.HlaAllele
import junit.framework.Assert.assertEquals
import org.junit.Test


class HlaComplexCoverageRankingTest {

    val victim = HlaComplexCoverageRanking(3, listOf(), listOf())

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
    fun testMostCommon() {
        val common = HlaComplexCoverage.create(listOf(a1, a3, b1, b3, c1, c2))
        val lessCommon = mutableListOf<HlaComplexCoverage>()

        lessCommon.add(HlaComplexCoverage.create(listOf(a1, a2, b1, b2, c1, c2)))
        lessCommon.add(HlaComplexCoverage.create(listOf(a1, a3, b1, b2, c1, c2)))
        lessCommon.add(HlaComplexCoverage.create(listOf(a1, a2, b1, b3, c1, c2)))
        lessCommon.add(HlaComplexCoverage.create(listOf(a1, a2, b1, b2, c1, c3)))

        for (lower in lessCommon) {
            val complexes = mutableListOf(common, lower).shuffled()
            val winner = HlaComplexCoverageRanking(3, listOf(a3.allele, b3.allele, c3.allele), listOf()).candidateRanking(complexes)[0]
            assertEquals(common, winner)
        }
    }

    @Test
    fun testLeastRecovered() {
        val common = listOf(a1.allele, a2.allele, a3.allele)
        val recovered = listOf(a1.allele)

        val victim1 = HlaComplexCoverage.create(listOf(a1, a2, b1, b2, c1, c2))
        val victim2 = HlaComplexCoverage.create(listOf(a1, a3, b1, b2, c1, c2))
        val victim3 = HlaComplexCoverage.create(listOf(a2, a3, b1, b2, c1, c2))
        val complexes = listOf(victim1, victim2, victim3).shuffled()

        val ranked = HlaComplexCoverageRanking(3, common, recovered).candidateRanking(complexes)
        val winner = ranked[0]
        assertEquals(victim3, winner)
    }

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
        val hetWithOneWild = HlaComplexCoverage.create(listOf(a1, a2, b1, b2, c1, c2.copy(uniqueCoverage = 9, wildCoverage = 1.0)))

        val hom = HlaComplexCoverage.create(listOf(a1, a1, b1, b1, c1, c1))
        val homWithOneWild = HlaComplexCoverage.create(listOf(a1, a1, b1, b1, c1, c1.copy(uniqueCoverage = 9, wildCoverage = 1.0)))
        val homWithTwoWild = HlaComplexCoverage.create(listOf(a1, a1, b1, b1.copy(uniqueCoverage = 9, wildCoverage = 1.0), c1, c1.copy(uniqueCoverage = 9, wildCoverage = 1.0)))

        assertEquals(hom, victim.candidateRanking(mutableListOf(het, hom).shuffled())[0])
        assertEquals(het, victim.candidateRanking(mutableListOf(het, homWithOneWild).shuffled())[0])
        assertEquals(homWithOneWild, victim.candidateRanking(mutableListOf(homWithOneWild, hetWithOneWild).shuffled())[0])
        assertEquals(hetWithOneWild, victim.candidateRanking(mutableListOf(homWithTwoWild, hetWithOneWild).shuffled())[0])
    }


}
