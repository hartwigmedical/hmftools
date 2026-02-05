package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.annotation.ImgtSequenceFile
import com.hartwig.hmftools.cider.genes.VJ
import com.hartwig.hmftools.common.bam.CigarUtils.cigarElementsFromStr
import com.hartwig.hmftools.common.genome.region.Strand
import kotlin.test.*

class SomaticHypermutationTest {
    @Test
    fun testDetermineFromVGeneIdentity()
    {
        assertEquals(SomaticHypermutationStatus.UNMUTATED, SomaticHypermutationStatus.determineFromVGeneIdentity(98.0))
        assertEquals(SomaticHypermutationStatus.UNMUTATED, SomaticHypermutationStatus.determineFromVGeneIdentity(98.01))
        assertEquals(SomaticHypermutationStatus.UNMUTATED, SomaticHypermutationStatus.determineFromVGeneIdentity(99.0))
        assertEquals(SomaticHypermutationStatus.UNMUTATED, SomaticHypermutationStatus.determineFromVGeneIdentity(99.99))
        assertEquals(SomaticHypermutationStatus.UNMUTATED, SomaticHypermutationStatus.determineFromVGeneIdentity(100.0))

        assertEquals(SomaticHypermutationStatus.MUTATED_BORDERLINE, SomaticHypermutationStatus.determineFromVGeneIdentity(97.0))
        assertEquals(SomaticHypermutationStatus.MUTATED_BORDERLINE, SomaticHypermutationStatus.determineFromVGeneIdentity(97.01))
        assertEquals(SomaticHypermutationStatus.MUTATED_BORDERLINE, SomaticHypermutationStatus.determineFromVGeneIdentity(97.99))

        assertEquals(SomaticHypermutationStatus.MUTATED, SomaticHypermutationStatus.determineFromVGeneIdentity(96.99))
        assertEquals(SomaticHypermutationStatus.MUTATED, SomaticHypermutationStatus.determineFromVGeneIdentity(96.0))
        assertEquals(SomaticHypermutationStatus.MUTATED, SomaticHypermutationStatus.determineFromVGeneIdentity(95.0))
        assertEquals(SomaticHypermutationStatus.MUTATED, SomaticHypermutationStatus.determineFromVGeneIdentity(90.0))
    }

    @Test
    fun testCompareVJRegionToImgtVForward()
    {
        //   sections:   |---V----||-anchor-||--CDR3--||-anchor-||---J----|
        //   query:        |-------------------------------------------|
        //   alignment:     |------------------|
        val layoutSeq = "GGGTTGGGGGCAAAAACAAACCCCCCCCCCAAAAAAAAAATTTTTTTTTT"
        //   IMGT             |------------------|
        //   ref:        -----                    -----
        //   alignment:     |------------------|
        val imgtSeq =   "AAAAAGGGGGAAAAAAAAAACCCCCAAAAA"
        //   compare:         |-------------|

        val anchorBoundary = 20
        val queryRange = 2 until 47
        val imgtRange = 5 until 25
        val alignment = Alignment(
            layoutSeq.substring(queryRange), 1 until 21,
            "test", 3 until 23, Strand.FORWARD,
            20 * 4, 25, cigarElementsFromStr("1S20M24S"),
            imgtSeq.length)
        val actual = compareVJRegionToImgt(
            layoutSeq,
            VJ.V,
            anchorBoundary,
            ImgtSequenceFile.Sequence("test", "1", imgtSeq, imgtRange),
            queryRange,
            alignment)
        val expected = ShmGeneComparison(15, 15, 100.0 * 13 / 15, 0)
        assertEquals(expected, actual)
    }

    @Test
    fun testCompareVJRegionToImgtVReverse()
    {
        // TODO
    }

    @Test
    fun testCompareVJRegionToImgtJForward()
    {
        // TODO
    }

    @Test
    fun testCompareVJRegionToImgtJReverse()
    {
        // TODO
    }
}