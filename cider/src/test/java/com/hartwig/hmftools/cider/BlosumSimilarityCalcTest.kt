package com.hartwig.hmftools.cider

import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertNotEquals

class BlosumSimilarityCalcTest
{
    @Test
    fun testFindSimilarityCalcSimple()
    {
        val templateDnaSeq = "CAACCTGAAGATTTTGCAACTTACTACTGT"
        val testDnaSeq = "CAACCTGAAGATGTTGGATTTTACTTCTGT"
        var score = BlosumSimilarityCalc.calcSimilarityScore(VJ.V, templateDnaSeq, testDnaSeq)
        assertEquals(2, score)

        // test that one is longer than the other, we should align them correctly.
        // for V align right, for J align left
        score = BlosumSimilarityCalc.calcSimilarityScore(VJ.V, "TGTTGT$templateDnaSeq", testDnaSeq)
        assertEquals(2, score)
        score = BlosumSimilarityCalc.calcSimilarityScore(VJ.J, "${templateDnaSeq}TGTTGT", testDnaSeq)
        assertEquals(2, score)
    }
}