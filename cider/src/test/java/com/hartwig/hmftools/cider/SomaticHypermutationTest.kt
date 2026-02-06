package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.annotation.ImgtSequenceFile
import com.hartwig.hmftools.cider.genes.VJ
import com.hartwig.hmftools.common.bam.CigarUtils.cigarElementsFromStr
import com.hartwig.hmftools.common.genome.region.Strand
import htsjdk.samtools.util.SequenceUtil.reverseComplement
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
        val queryAlignRange = 1 until 21
        val refAlignRange = 3 until 23
        val alignment = Alignment(
            layoutSeq.substring(queryRange), queryAlignRange,
            "test", refAlignRange, Strand.FORWARD,
            0, 27, cigarElementsFromStr("1S20M24S"),
            imgtSeq.length)
        val actual = compareVJRegionToImgt(
            layoutSeq,
            VJ.V,
            anchorBoundary,
            ImgtSequenceFile.Sequence("test", "1", imgtSeq, 5, 5),
            queryRange,
            alignment)
        val expected = ShmGeneComparison(15, 15, 100.0 * 13 / 15, 0)
        assertEquals(expected, actual)
    }

    @Test
    fun testCompareVJRegionToImgtVReverse()
    {
        //   sections:   |---V----||-anchor-||--CDR3--||-anchor-||---J----|
        //   query:        |-------------------------------------------|
        //   alignment:     |------------------|
        val layoutSeq = "GGGTTGGGGGCAAAAACAAACCCCCCCCCCAAAAAAAAAATTTTTTTTTT"
        //   IMGT             |------------------|
        //   ref:        -----                    -----
        //   alignment:     |------------------|
        var imgtSeq =   "AAAAAGGGGGAAAAAAAAAACCCCCAAAAA"
        //   compare:         |-------------|

        imgtSeq = reverseComplement(imgtSeq)

        val anchorBoundary = 20
        val queryRange = 2 until 47
        val queryAlignRange = 1 until 21
        val refAlignRange = imgtSeq.length - 23 until imgtSeq.length - 3
        val alignment = Alignment(
            layoutSeq.substring(queryRange), queryAlignRange,
            "test", refAlignRange, Strand.REVERSE,
            0, 27, cigarElementsFromStr("24S20M1S"),
            imgtSeq.length)
        val actual = compareVJRegionToImgt(
            layoutSeq,
            VJ.V,
            anchorBoundary,
            ImgtSequenceFile.Sequence("test", "1", imgtSeq, 5, 5),
            queryRange,
            alignment)
        val expected = ShmGeneComparison(15, 15, 100.0 * 13 / 15, 0)
        assertEquals(expected, actual)
    }

    @Test
    fun testCompareVJRegionToImgtJForward()
    {
        //   sections:   |---V----||-anchor-||--CDR3--||-anchor-||---J----|
        //   query:        |--------------------------------------------|
        //   alignment:                             |------------------|
        val layoutSeq = "GGGGGGGGGGAAAAAAAAAACCCCCCCTCTCAAAAAAAAGTTTTTTTTTT"
        //   IMGT                                 |------------------|
        //   ref:                            -----                    -----
        //   alignment:                             |------------------|
        val imgtSeq =                       "AAAAACCCCCAAAAAAAAAATTTTTAAAAA"
        //   compare:                                  |-------------|

        val anchorBoundary = 30
        val queryRange = 2 until 48
        val queryAlignRange = 25 until 45
        val refAlignRange = 7 until 27
        val alignment = Alignment(
            layoutSeq.substring(queryRange), queryAlignRange,
            "test", refAlignRange, Strand.FORWARD,
            0, 27, cigarElementsFromStr("25S20M1S"),
            imgtSeq.length)
        val actual = compareVJRegionToImgt(
            layoutSeq,
            VJ.J,
            anchorBoundary,
            ImgtSequenceFile.Sequence("test", "1", imgtSeq, 5, 5),
            queryRange,
            alignment)
        val expected = ShmGeneComparison(15, 15, 100.0 * 13 / 15, 0)
        assertEquals(expected, actual)
    }

    @Test
    fun testCompareVJRegionToImgtJReverse()
    {
        //   sections:   |---V----||-anchor-||--CDR3--||-anchor-||---J----|
        //   query:        |--------------------------------------------|
        //   alignment:                             |------------------|
        val layoutSeq = "GGGGGGGGGGAAAAAAAAAACCCCCCCTCTCAAAAAAAAGTTTTTTTTTT"
        //   IMGT                                 |------------------|
        //   ref:                            -----                    -----
        //   alignment:                             |------------------|
        var imgtSeq =                       "AAAAACCCCCAAAAAAAAAATTTTTAAAAA"
        //   compare:                                  |-------------|

        imgtSeq = reverseComplement(imgtSeq)

        val anchorBoundary = 30
        val queryRange = 2 until 48
        val queryAlignRange = 25 until 45
        val refAlignRange = imgtSeq.length - 27 until imgtSeq.length - 7
        val alignment = Alignment(
            layoutSeq.substring(queryRange), queryAlignRange,
            "test", refAlignRange, Strand.REVERSE,
            0, 27, cigarElementsFromStr("1S20M25S"),
            imgtSeq.length)
        val actual = compareVJRegionToImgt(
            layoutSeq,
            VJ.J,
            anchorBoundary,
            ImgtSequenceFile.Sequence("test", "1", imgtSeq, 5, 5),
            queryRange,
            alignment)
        val expected = ShmGeneComparison(15, 15, 100.0 * 13 / 15, 0)
        assertEquals(expected, actual)
    }
}