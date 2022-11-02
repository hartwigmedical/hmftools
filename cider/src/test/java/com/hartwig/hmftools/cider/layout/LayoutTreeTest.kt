/*

package com.hartwig.hmftools.cider.layout

import com.hartwig.hmftools.cider.ReadKey
import htsjdk.samtools.SAMUtils
import org.apache.logging.log4j.Level
import org.junit.Before
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class LayoutTreeTest
{
    @Before
    fun setUp()
    {
        org.apache.logging.log4j.core.config.Configurator.setRootLevel(Level.TRACE);
    }

    @Test
    fun testAddReadSimple()
    {
        // test very simple add high qual reads
        val layoutTree = LayoutTree(MIN_BASE_QUALITY, MIN_OVERLAP_BASES)
        val seq1 = "CAGGT"
        val baseQual1 = SAMUtils.fastqToPhred("FFFFF") // F is 37, : is 25

        val read1 = LayoutTree.Read("read1", ReadKey("read1", true), seq1, baseQual1, 0)

        // add it to the layout tree
        assertTrue(layoutTree.tryAddRead(read1))

        // now should have 5 levls
        assertEquals(5, layoutTree.numLevels)

        // check we can get the same layout
        var layouts = layoutTree.buildReadLayouts()

        assertEquals(1, layouts.size)
        assertEquals(seq1, layouts.first().consensusSequence())
        assertEquals("11111", layouts.first().highQualSupportString())

        // add another read of same thing, but has 1 low qual base
        val baseQual2 = SAMUtils.fastqToPhred("FF:FF") // F is 37, : is 25
        val read2 = LayoutTree.Read("read2", ReadKey("read2", true), seq1, baseQual2, 0)

        // add it to the layout tree
        assertTrue(layoutTree.tryAddRead(read2))
        layouts = layoutTree.buildReadLayouts()
        assertEquals(1, layouts.size)
        assertEquals(seq1, layouts.first().consensusSequence())
        assertEquals("22122", layouts.first().highQualSupportString())
    }

    @Test
    fun testAddReadSimpleBranch()
    {
        // test very simple add high qual reads
        val layoutTree = LayoutTree(MIN_BASE_QUALITY, MIN_OVERLAP_BASES)
        val seq1 = "CAGGT"
        val baseQual1 = SAMUtils.fastqToPhred("FFFFF") // F is 37, : is 25

        val read1 = LayoutTree.Read("read1", ReadKey("read1", true), seq1, baseQual1, 0)

        // seq 2 is different for the middle letter
        val seq2 = "CATGT"
        val baseQual2 = SAMUtils.fastqToPhred("FFFFF") // F is 37, : is 25
        val read2 = LayoutTree.Read("read2", ReadKey("read2", true), seq2, baseQual2, 0)

        // add it to the layout tree
        assertTrue(layoutTree.tryAddRead(read1))
        assertTrue(layoutTree.tryAddRead(read2))
        val layouts = layoutTree.buildReadLayouts().sortedBy({ l -> l.consensusSequence() })
        assertEquals(2, layouts.size)
        assertEquals(seq1, layouts[0].consensusSequence())
        assertEquals("11111", layouts[0].highQualSupportString())
        assertEquals(seq2, layouts[1].consensusSequence())
        assertEquals("11111", layouts[1].highQualSupportString())
    }

    @Test
    fun testAddReadLowQualMismatch()
    {
        // if we encounter a low qual mismatch, and can put it into other branch it is accepted

        val layoutTree = LayoutTree(MIN_BASE_QUALITY, MIN_OVERLAP_BASES)
        val seq1 = "CAGGT"
        val baseQual1 = SAMUtils.fastqToPhred("FFFFF") // F is 37, : is 25

        val read1 = LayoutTree.Read("read1", ReadKey("read1", true), seq1, baseQual1, 0)

        // seq 2 is different for the middle letter
        val seq2 = "CATAT"
        val baseQual2 = SAMUtils.fastqToPhred("FF::F") // F is 37, : is 25
        val read2 = LayoutTree.Read("read2", ReadKey("read2", true), seq2, baseQual2, 0)

        // add it to the layout tree
        assertTrue(layoutTree.tryAddRead(read1))
        assertTrue(layoutTree.tryAddRead(read2))

        // since both match the high quality one we only get one layout
        val layouts = layoutTree.buildReadLayouts().sortedBy({ l -> l.consensusSequence() })
        assertEquals(1, layouts.size)
        assertEquals(seq1, layouts[0].consensusSequence())
        assertEquals("22112", layouts[0].highQualSupportString())
    }

    @Test
    fun testAddReadLowQualMismatch1()
    {
        // we first add the low qual one and then high qual one

        val layoutTree = LayoutTree(MIN_BASE_QUALITY, MIN_OVERLAP_BASES)
        val seq1 = "CATAT"
        val baseQual1 = SAMUtils.fastqToPhred("FF::F") // F is 37, : is 25

        val read1 = LayoutTree.Read("read1", ReadKey("read1", true), seq1, baseQual1, 0)

        // seq 2 is different for the middle letter
        val seq2 = "CAGCT"
        val baseQual2 = SAMUtils.fastqToPhred("FFFFF") // F is 37, : is 25
        val read2 = LayoutTree.Read("read2", ReadKey("read2", true), seq2, baseQual2, 0)

        // add it to the layout tree
        assertTrue(layoutTree.tryAddRead(read1))
        assertTrue(layoutTree.tryAddRead(read2))

        // since both match the high quality one we only get one layout
        val layouts = layoutTree.buildReadLayouts().sortedBy({ l -> l.consensusSequence() })
        assertEquals(1, layouts.size)
        assertEquals(seq2, layouts[0].consensusSequence())
        assertEquals("22112", layouts[0].highQualSupportString())
    }

    @Test
    fun testLowQualBranch()
    {
        // we encounter a low qual branch that needs to be confirmed

        val layoutTree = LayoutTree(MIN_BASE_QUALITY, MIN_OVERLAP_BASES)

        val seq1 = "CAGCT"
        val baseQual1 = SAMUtils.fastqToPhred("FFFFF") // F is 37, : is 25
        val read1 = LayoutTree.Read("read1", ReadKey("read1", true), seq1, baseQual1, 0)

        // seq 2 has a low qual diff in 3rd and high qual diff at 5th
        val seq2 = "CATCA"
        val baseQual2 = SAMUtils.fastqToPhred("FF:FF") // F is 37, : is 25
        val read2 = LayoutTree.Read("read2", ReadKey("read2", true), seq2, baseQual2, 0)

        // seq 3 is tricky, it actually flips the 3rd base to G, so the branch needs to be moved
        val seq3 = "CAGCA"
        val baseQual3 = SAMUtils.fastqToPhred("FFFFF") // F is 37, : is 25
        val read3 = LayoutTree.Read("read3", ReadKey("read3", true), seq3, baseQual3, 0)

        // add to the layout tree
        assertTrue(layoutTree.tryAddRead(read1))
        assertTrue(layoutTree.tryAddRead(read2))
        assertTrue(layoutTree.tryAddRead(read3))

        // since both match the high quality one we only get one layout
        val layouts = layoutTree.buildReadLayouts().sortedBy({ l -> l.consensusSequence() })
        assertEquals(2, layouts.size)
        assertEquals(seq3, layouts[0].consensusSequence())
        assertEquals("22122", layouts[0].highQualSupportString())
        assertEquals(seq1, layouts[1].consensusSequence())
        assertEquals("11111", layouts[1].highQualSupportString())
    }

    @Test
    fun testLowQualBranchExtend()
    {
        // we encounter a low qual branch that needs to be confirmed

        val layoutTree = LayoutTree(MIN_BASE_QUALITY, MIN_OVERLAP_BASES)

        val seq1 = "CAGCT"
        val baseQual1 = SAMUtils.fastqToPhred("FF:FF") // F is 37, : is 25
        val read1 = LayoutTree.Read("read1", ReadKey("read1", true), seq1, baseQual1, 0)

        // seq 2 has a low qual diff in 3rd and high qual diff at 5th
        val seq2 = "CATCA"
        val baseQual2 = SAMUtils.fastqToPhred("FF:FF") // F is 37, : is 25
        val read2 = LayoutTree.Read("read2", ReadKey("read2", true), seq2, baseQual2, 0)

        // this one actually can merge with seq2, and extend it such that end result is CAGCAT
        val seq3 = "CAGTAT"
        val baseQual3 = SAMUtils.fastqToPhred("FFF:FF") // F is 37, : is 25
        val read3 = LayoutTree.Read("read3", ReadKey("read3", true), seq3, baseQual3, 0)

        // add to the layout tree
        assertTrue(layoutTree.tryAddRead(read1))
        assertTrue(layoutTree.tryAddRead(read2))
        assertTrue(layoutTree.tryAddRead(read3))

        // since both match the high quality one we only get one layout
        val layouts = layoutTree.buildReadLayouts().sortedBy({ l -> l.consensusSequence() })
        assertEquals(2, layouts.size)
        assertEquals("CAGCAT", layouts[0].consensusSequence())
        assertEquals("221121", layouts[0].highQualSupportString())
        assertEquals(seq1, layouts[1].consensusSequence())
        assertEquals("11011", layouts[1].highQualSupportString())
    }

    companion object
    {
        const val MIN_BASE_QUALITY = 30.toByte()
        const val MIN_OVERLAP_BASES = 3
    }
}

 */