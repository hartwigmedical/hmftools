package com.hartwig.hmftools.paddle.dnds

import com.hartwig.hmftools.paddle.Impact
import org.junit.Assert
import org.junit.Assert.assertEquals
import org.junit.Test

class DndsMutationTest {

    companion object {

        const val GENE = "AR"

        fun dndsMutation(gene: String, hotspot: Boolean, biallelic: Boolean, repeatCount: Int, impact: Impact): DndsMutation {
            val biallelicString = biallelic.toString()
            val hotspotString = hotspot.toString()
            val dndImpactString = when (impact) {
                Impact.INFRAME, Impact.FRAMESHIFT -> "no-SNV"
                Impact.MISSENSE -> "Missense"
                Impact.SPLICE -> "Essential_Splice"
                Impact.NONSENSE -> "Nonsense"
                Impact.SYNONYMOUS -> "Synonymous"
                Impact.UNKNOWN -> "NA"
            }

            val canonicalCodingEffectString = when (impact) {
                Impact.INFRAME -> "MISSENSE"
                Impact.FRAMESHIFT -> "NONSENSE_OR_FRAMESHIFT"
                else -> "NA"
            }

            val result = DndsMutation.fromString("SAMPLE\tX\t66766353\tT\tTGGTGGCGGC\t${canonicalCodingEffectString}\t${canonicalCodingEffectString}\t${repeatCount}\t${biallelicString}\t${hotspotString}\t${gene}\t${dndImpactString}")
            assertEquals(impact, result.impact)
            return result
        }
    }

    private val inframe = DndsMutation.fromString("SAMPLE\tX\t66766353\tT\tTGGTGGCGGC\tMISSENSE\tMISSENSE\t2\tfalse\tfalse\tAR\tno-SNV")
    private val worstInframe = DndsMutation.fromString("SAMPLE\t19\t45853980\tGG\tAA\tMISSENSE\tNONE\t2\tfalse\tfalse\tKLC3\tno-SNV")
    private val canonicalInframe = DndsMutation.fromString("SAMPLE\t8\t145743003\tGCCGCCTCCACGTCGT\tG\tNONE\tMISSENSE\t2\tfalse\tfalse\tRECQL4\tno-SNV")

    private val frameshift = DndsMutation.fromString("SAMPLE\tX\t66766423\tGC\tG\tNONSENSE_OR_FRAMESHIFT\tNONSENSE_OR_FRAMESHIFT\t5\tfalse\tfalse\tAR\tno-SNV")
    private val worstFrameshift = DndsMutation.fromString("SAMPLE\tX\t70357238\tAG\tA\tNONSENSE_OR_FRAMESHIFT\tNONE\t0\tfalse\tfalse\tMED12\tno-SNV")

    private val missense = DndsMutation.fromString("SAMPLE\tX\t66766561\tG\tT\tMISSENSE\tMISSENSE\t3\tfalse\tfalse\tAR\tMissense")
    private val nonsense = DndsMutation.fromString("SAMPLE\tX\t66765547\tG\tT\tNONSENSE_OR_FRAMESHIFT\tNONSENSE_OR_FRAMESHIFT\t0\tfalse\tfalse\tAR\tNonsense")
    private val splice = DndsMutation.fromString("SAMPLE\t19\t7795746\tC\tA\tSPLICE\tSPLICE\t2\tfalse\tfalse\tCLEC4G\tEssential_Splice")
    private val synonymous = DndsMutation.fromString("SAMPLE\tX\t66766206\tC\tT\tSYNONYMOUS\tSYNONYMOUS\t2\tfalse\tfalse\tAR\tSynonymous")
    private val unknown = DndsMutation.fromString("SAMPLE\tX\t53576468\tTG\tT\tNONE\tNONE\t2\tfalse\tfalse\tHUWE1\tno-SNV")

    private val hotspotFrameshift = DndsMutation.fromString("SAMPLE\tX\t66766423\tGC\tG\tNONSENSE_OR_FRAMESHIFT\tNONSENSE_OR_FRAMESHIFT\t5\tfalse\ttrue\tAR\tno-SNV")
    private val biallelicNonsense = DndsMutation.fromString("SAMPLE\tX\t66765547\tG\tT\tNONSENSE_OR_FRAMESHIFT\tNONSENSE_OR_FRAMESHIFT\t0\ttrue\tfalse\tAR\tNonsense")

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
        list.sortWith(DndsMutationComparator { x -> x.isKnownOncoDriver })

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
        list.sortWith(DndsMutationComparator { x -> x.isKnownTsgDriver })

        assertEquals(Impact.FRAMESHIFT, list[0].impact)
        assertEquals(Impact.NONSENSE, list[1].impact)
        assertEquals(Impact.INFRAME, list[2].impact)
    }

    @Test
    fun testSortByHotspot() {
        val nonsenseHotspot = dndsMutation(GENE, true, false, 0, Impact.NONSENSE)
        val inframe = dndsMutation(GENE, false, false, 0, Impact.INFRAME)
        val list = mutableListOf(inframe, nonsenseHotspot)
        list.shuffle()
        list.sortWith(DndsMutationComparator { x -> x.isKnownOncoDriver })
        assertEquals(Impact.NONSENSE, list[0].impact)
        assertEquals(Impact.INFRAME, list[1].impact)
    }

    @Test
    fun testIsHotspot() {
        Assert.assertFalse(dndsMutation(GENE, true, false, 0, Impact.SYNONYMOUS).isHotspot)
        Assert.assertFalse(dndsMutation(GENE, true, false, 0, Impact.UNKNOWN).isHotspot)

        Assert.assertTrue(dndsMutation(GENE, true, false, 0, Impact.MISSENSE).isHotspot)
        Assert.assertTrue(dndsMutation(GENE, true, false, 0, Impact.NONSENSE).isHotspot)
        Assert.assertTrue(dndsMutation(GENE, true, false, 0, Impact.SPLICE).isHotspot)
        Assert.assertTrue(dndsMutation(GENE, true, false, 0, Impact.INFRAME).isHotspot)
        Assert.assertTrue(dndsMutation(GENE, true, false, 0, Impact.FRAMESHIFT).isHotspot)
    }

    @Test
    fun testIsBiallelic() {
        Assert.assertFalse(dndsMutation(GENE, false, true, 0, Impact.MISSENSE).isBiallelic)
        Assert.assertFalse(dndsMutation(GENE, false, true, 0, Impact.SYNONYMOUS).isBiallelic)
        Assert.assertFalse(dndsMutation(GENE, false, true, 0, Impact.UNKNOWN).isBiallelic)

        Assert.assertTrue(dndsMutation(GENE, false, true, 0, Impact.NONSENSE).isBiallelic)
        Assert.assertTrue(dndsMutation(GENE, false, true, 0, Impact.SPLICE).isBiallelic)
        Assert.assertTrue(dndsMutation(GENE, false, true, 0, Impact.INFRAME).isBiallelic)
        Assert.assertTrue(dndsMutation(GENE, false, true, 0, Impact.FRAMESHIFT).isBiallelic)
    }

    @Test
    fun testKnownTsg() {
        Assert.assertFalse(dndsMutation(GENE, false, false, 0, Impact.INFRAME).isKnownTsgDriver)
        Assert.assertTrue(dndsMutation(GENE, true, false, 0, Impact.INFRAME).isKnownTsgDriver)
        Assert.assertTrue(dndsMutation(GENE, true, true, 0, Impact.INFRAME).isKnownTsgDriver)
    }

    @Test
    fun testKnownOnco() {
        Assert.assertFalse(dndsMutation(GENE, false, true, 8, Impact.INFRAME).isKnownOncoDriver)
        Assert.assertTrue(dndsMutation(GENE, true, false, 8, Impact.INFRAME).isKnownOncoDriver)
        Assert.assertTrue(dndsMutation(GENE, false, false, 7, Impact.INFRAME).isKnownOncoDriver)
        Assert.assertFalse(dndsMutation(GENE, false, false, 7, Impact.FRAMESHIFT).isKnownOncoDriver)
    }
}