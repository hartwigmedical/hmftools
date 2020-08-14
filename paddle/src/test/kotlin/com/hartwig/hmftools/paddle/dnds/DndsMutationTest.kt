package com.hartwig.hmftools.paddle.dnds

import junit.framework.Assert.assertEquals
import org.junit.Test

class DndsMutationTest {

    private val inframe = DndsMutation.fromString("SAMPLE\tX\t66766353\tT\tTGGTGGCGGC\tMISSENSE\tMISSENSE\t2\t0\tNON_HOTSPOT\tAR\tno-SNV")
    private val worstInframe = DndsMutation.fromString("SAMPLE\t19\t45853980\tGG\tAA\tMISSENSE\tNONE\t2\t0\tNON_HOTSPOT\tKLC3\tno-SNV")
    private val canonicalInframe = DndsMutation.fromString("SAMPLE\t8\t145743003\tGCCGCCTCCACGTCGT\tG\tNONE\tMISSENSE\t2\t0\tNON_HOTSPOT\tRECQL4\tno-SNV")

    private val frameshift = DndsMutation.fromString("SAMPLE\tX\t66766423\tGC\tG\tNONSENSE_OR_FRAMESHIFT\tNONSENSE_OR_FRAMESHIFT\t5\t0\tNON_HOTSPOT\tAR\tno-SNV")
    private val worstFrameshift = DndsMutation.fromString("SAMPLE\tX\t70357238\tAG\tA\tNONSENSE_OR_FRAMESHIFT\tNONE\t0\t0\tNON_HOTSPOT\tMED12\tno-SNV")

    private val missense = DndsMutation.fromString("SAMPLE\tX\t66766561\tG\tT\tMISSENSE\tMISSENSE\t3\t0\tNON_HOTSPOT\tAR\tMissense")
    private val nonsense = DndsMutation.fromString("SAMPLE\tX\t66765547\tG\tT\tNONSENSE_OR_FRAMESHIFT\tNONSENSE_OR_FRAMESHIFT\t0\t0\tNON_HOTSPOT\tAR\tNonsense")
    private val splice = DndsMutation.fromString("SAMPLE\t19\t7795746\tC\tA\tSPLICE\tSPLICE\t2\t0\tNON_HOTSPOT\tCLEC4G\tEssential_Splice")
    private val synonymous = DndsMutation.fromString("SAMPLE\tX\t66766206\tC\tT\tSYNONYMOUS\tSYNONYMOUS\t2\t0\tNON_HOTSPOT\tAR\tSynonymous")
    private val unknown = DndsMutation.fromString("SAMPLE\tX\t53576468\tTG\tT\tNONE\tNONE\t2\t0\tNON_HOTSPOT\tHUWE1\tno-SNV")

    private val hotspotFrameshift = DndsMutation.fromString("SAMPLE\tX\t66766423\tGC\tG\tNONSENSE_OR_FRAMESHIFT\tNONSENSE_OR_FRAMESHIFT\t5\t0\tHOTSPOT\tAR\tno-SNV")
    private val biallelicNonsense = DndsMutation.fromString("SAMPLE\tX\t66765547\tG\tT\tNONSENSE_OR_FRAMESHIFT\tNONSENSE_OR_FRAMESHIFT\t0\t1\tNON_HOTSPOT\tAR\tNonsense")

    @Test
    fun testImpact() {
        assertEquals(Impact.MISSENSE, missense.impact)
        assertEquals(Impact.FRAMESHIFT, frameshift.impact)
        assertEquals(Impact.FRAMESHIFT, worstFrameshift.impact)
        assertEquals(Impact.NONSENSE, nonsense.impact)
        assertEquals(Impact.SPLICE, splice.impact)
        assertEquals(Impact.SYNONYMOUS, synonymous.impact)
        assertEquals(Impact.INFRAME, inframe.impact)
        assertEquals(Impact.INFRAME, canonicalInframe.impact)
        assertEquals(Impact.INFRAME, worstInframe.impact)
        assertEquals(Impact.UNKNOWN, unknown.impact)
    }

    @Test
    fun testComparator() {
        val list = mutableListOf(
                missense, frameshift, inframe, nonsense, splice, synonymous, worstInframe, worstFrameshift, canonicalInframe, unknown)
        list.shuffle()
        list.sortWith(DndsMutationComparator(false))

        assertEquals(Impact.INFRAME, list[0].impact)
        assertEquals(Impact.INFRAME, list[1].impact)
        assertEquals(Impact.INFRAME, list[2].impact)
        assertEquals(Impact.MISSENSE, list[3].impact)
        assertEquals(Impact.SPLICE, list[4].impact)
        assertEquals(Impact.NONSENSE, list[5].impact)
        assertEquals(Impact.FRAMESHIFT, list[6].impact)
        assertEquals(Impact.FRAMESHIFT, list[7].impact)
        assertEquals(Impact.SYNONYMOUS, list[8].impact)
        assertEquals(Impact.UNKNOWN, list[9].impact)
    }

    @Test
    fun testHotspotBiallelicSorting() {
        val list = mutableListOf(inframe, hotspotFrameshift, biallelicNonsense)
        list.shuffle()
        list.sortWith(DndsMutationComparator(true))

        assertEquals(Impact.FRAMESHIFT, list[0].impact)
        assertEquals(Impact.NONSENSE, list[1].impact)
        assertEquals(Impact.INFRAME, list[2].impact)

    }



}