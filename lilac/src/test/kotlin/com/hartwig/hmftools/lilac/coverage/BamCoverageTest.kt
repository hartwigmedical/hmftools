package com.hartwig.hmftools.lilac.coverage

import com.hartwig.hmftools.lilac.dna.aminoAcids
import com.hartwig.hmftools.lilac.hla.HlaAllele
import org.junit.Assert
import org.junit.Test

class BamCoverageTest {

    @Test
    fun testProcessDna() {

        val randomDna = "CATCATCATCAT"
        val commonDna = "AATATAGCCATCATTATGGAAAAATTAAACCTAACAGTATTA"
        val uniqueDNA1 = "CGTAACAGTTGGATTGGGCAGGATAGCATTGTTTTAACCCGG"
        val uniqueDNA2 = "AGGTGTGGGCTGGGCTGGTCAGGTGTGGGACGGGCTGGTCAG"

        val commonAA = commonDna.aminoAcids()
        val uniqueAA1 = uniqueDNA1.aminoAcids()
        val uniqueAA2 = uniqueDNA2.aminoAcids()

        val proteinCoverage1 = ProteinCoverage(HlaAllele("A*01:01"), listOf(commonAA, uniqueAA1))
        val proteinCoverage2 = ProteinCoverage(HlaAllele("B*01:02"), listOf(commonAA, uniqueAA2))
        val victim = BamCoverage(3, setOf(proteinCoverage1, proteinCoverage2))

        assertStuff(intArrayOf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), proteinCoverage1.map[commonAA]!!)
        assertStuff(intArrayOf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), proteinCoverage1.map[uniqueAA1]!!)
        assertStuff(intArrayOf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), proteinCoverage2.map[commonAA]!!)
        assertStuff(intArrayOf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), proteinCoverage2.map[uniqueAA2]!!)

        // Match isn't big enough
        victim.processDna(randomDna + uniqueDNA1.startCodons(2))
        assertStuff(intArrayOf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), proteinCoverage1.map[commonAA]!!)
        assertStuff(intArrayOf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), proteinCoverage1.map[uniqueAA1]!!)
        assertStuff(intArrayOf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), proteinCoverage2.map[commonAA]!!)
        assertStuff(intArrayOf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), proteinCoverage2.map[uniqueAA2]!!)


        // Match start of unique 1
        victim.processDna(randomDna + uniqueDNA1.startCodons(3))
        assertStuff(intArrayOf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), proteinCoverage1.map[commonAA]!!)
        assertStuff(intArrayOf(1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), proteinCoverage1.map[uniqueAA1]!!)
        assertStuff(intArrayOf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), proteinCoverage2.map[commonAA]!!)
        assertStuff(intArrayOf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), proteinCoverage2.map[uniqueAA2]!!)


        // Match end of unique 2
        victim.processDna(uniqueDNA2.endCodons(4) + randomDna)
        assertStuff(intArrayOf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), proteinCoverage1.map[commonAA]!!)
        assertStuff(intArrayOf(1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), proteinCoverage1.map[uniqueAA1]!!)
        assertStuff(intArrayOf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), proteinCoverage2.map[commonAA]!!)
        assertStuff(intArrayOf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1), proteinCoverage2.map[uniqueAA2]!!)

        // Match start of common
        victim.processDna(randomDna + commonDna.startCodons(6))
        assertStuff(intArrayOf(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0), proteinCoverage1.map[commonAA]!!)
        assertStuff(intArrayOf(1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), proteinCoverage1.map[uniqueAA1]!!)
        assertStuff(intArrayOf(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0), proteinCoverage2.map[commonAA]!!)
        assertStuff(intArrayOf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1), proteinCoverage2.map[uniqueAA2]!!)

        // Match end of common
        victim.processDna(commonDna.endCodons(10) + randomDna)
        assertStuff(intArrayOf(1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1), proteinCoverage1.map[commonAA]!!)
        assertStuff(intArrayOf(1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), proteinCoverage1.map[uniqueAA1]!!)
        assertStuff(intArrayOf(1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1), proteinCoverage2.map[commonAA]!!)
        assertStuff(intArrayOf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1), proteinCoverage2.map[uniqueAA2]!!)

        // Match all of common
        victim.processDna(randomDna + commonDna + randomDna)
        assertStuff(intArrayOf(2, 2, 2, 2, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2), proteinCoverage1.map[commonAA]!!)
        assertStuff(intArrayOf(1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), proteinCoverage1.map[uniqueAA1]!!)
        assertStuff(intArrayOf(2, 2, 2, 2, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2), proteinCoverage2.map[commonAA]!!)
        assertStuff(intArrayOf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1), proteinCoverage2.map[uniqueAA2]!!)

    }

    private fun assertStuff(expected: IntArray, exon: ExonCoverage) {
        val expectedString = expected.joinToString { it.toString() }
        val actualString = exon.coverage.joinToString { it.toString() }
        Assert.assertEquals(expectedString, actualString)
    }

    private fun String.endCodons(codons: Int): String {
        return this.substring(this.length - codons * 3, this.length)
    }

    private fun String.startCodons(codons: Int): String {
        return  this.substring(0, codons * 3 + 1)
    }

}