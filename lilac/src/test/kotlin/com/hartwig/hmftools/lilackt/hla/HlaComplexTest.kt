package com.hartwig.hmftools.lilackt.hla

import com.hartwig.hmftools.lilackt.coverage.HlaComplex
import junit.framework.Assert.assertEquals
import org.junit.Test
import java.util.*

class HlaComplexTest {

    val a11 = createA(1,1)
    val a12 = createA(1,2)
    val a13 = createA(1,3)
    val a21 = createA(2,1)
    val a22 = createA(2,2)
    val a23 = createA(2,3)
    val a31 = createA(3,1)

    val all = listOf(a11, a12, a13, a21, a22, a23, a31)

    val b11 = createB(1,1)
    val b12 = createB(1,2)
    val b13 = createB(1,3)

    @Test
    fun testTwoConfirmedProtein() {
        val confirmedProtein = listOf(a11, a12)
        val confirmedGroup = confirmedProtein.map { it.asAlleleGroup() }.distinct()

        val result = HlaComplex.gene("A", confirmedGroup, confirmedProtein, all)
        assertEquals(confirmedProtein, result[0].alleles)
    }

    @Test
    fun testOneConfirmedGroup() {
        val confirmedGroup = listOf(a11.asAlleleGroup())

        val result = HlaComplex.gene("A", confirmedGroup, Collections.emptyList(), all)
        for (hlaComplex in result) {
            println(hlaComplex)
        }

    }

    @Test
    fun testNoConfirmed() {
        val result = HlaComplex.gene("A", Collections.emptyList(), Collections.emptyList(), all)
        for (hlaComplex in result) {
            println(hlaComplex)
        }

    }

    @Test
    fun testCombineComplexes() {
        val complexA1 = HlaComplex(listOf(a11, a12))
        val complexA2 = HlaComplex(listOf(a11, a13))
        val complexA3 = HlaComplex(listOf(a12, a13))
        val complexA = listOf(complexA1, complexA2, complexA3)

        val complexB1 = HlaComplex(listOf(b11, b12))
        val complexB2 = HlaComplex(listOf(b11, b13))
        val complexB = listOf(complexB1, complexB2)

        val result = HlaComplex.combineComplexes(complexA, complexB)
        for (hlaComplex in result) {
            println(hlaComplex.alleles)
        }

    }

    private fun createA(group: Int, protein: Int): HlaAllele {
        return HlaAllele("A", group.toString(), protein.toString(), "", "")
    }

    private fun createB(group: Int, protein: Int): HlaAllele {
        return HlaAllele("B", group.toString(), protein.toString(), "", "")
    }

}