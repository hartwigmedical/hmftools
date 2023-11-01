package com.hartwig.hmftools.cider.layout

import com.hartwig.hmftools.cider.ReadKey
import htsjdk.samtools.SAMUtils
import org.junit.Before
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class LayoutForestTest
{
    @Before
    fun setUp()
    {
        // org.apache.logging.log4j.core.config.Configurator.setRootLevel(org.apache.logging.log4j.Level.TRACE)
    }

    val layoutReadCreateFunc = { read: LayoutForest.Read -> TestLayoutRead(
        read.source as String , ReadKey(read.source as String, true),
        read.sequence, read.baseQualities, read.alignedPosition) }

    @Test
    fun testAddReadSimple1()
    {
        // test very simple add high qual reads
        val layoutForest = LayoutForest(MIN_BASE_QUALITY, MIN_OVERLAP_BASES, MIN_SUPPORT_TO_SEAL_NODE)
        val seq1 = "CAGGT"
        val baseQual1 = SAMUtils.fastqToPhred("FFFFF") // F is 37, : is 25

        val read1 = LayoutForest.Read("read1", seq1.toByteArray(), baseQual1, 0)

        // add it to the layout tree
        assertTrue(layoutForest.tryAddRead(read1))

        // now should have 5 levls
        assertEquals(5, layoutForest.numLevels)

        // check we can get the same layout
        var layouts = layoutForest.buildReadLayouts(layoutReadCreateFunc)

        assertEquals(1, layouts.size)
        assertEquals(seq1, layouts.first().consensusSequenceString())
        assertEquals("11111", layouts.first().highQualSupportString())

        // add another read of same thing, but has 1 low qual base
        val baseQual2 = SAMUtils.fastqToPhred("FF:FF") // F is 37, : is 25
        val read2 = LayoutForest.Read("read2", seq1.toByteArray(), baseQual2, 0)

        // add it to the layout tree
        assertTrue(layoutForest.tryAddRead(read2))
        layouts = layoutForest.buildReadLayouts(layoutReadCreateFunc)
        assertEquals(1, layouts.size)
        assertEquals(seq1, layouts.first().consensusSequenceString())
        assertEquals("22122", layouts.first().highQualSupportString())
    }

    @Test
    fun testAddReadWrongOrder()
    {
        // test very simple add high qual reads
        val layoutForest = LayoutForest(MIN_BASE_QUALITY, MIN_OVERLAP_BASES, MIN_SUPPORT_TO_SEAL_NODE)

        layoutForest.tryAddRead(LayoutForest.Read("r2", 
            "CAGCT".toByteArray(), SAMUtils.fastqToPhred("FFFFF"), 4))

        layoutForest.tryAddRead(LayoutForest.Read("r1", 
            "GCAGCT".toByteArray(), SAMUtils.fastqToPhred("FFFFFF"), 5))

        // now should have 6 levls
        assertEquals(6, layoutForest.numLevels)

        // check we can get the same layout
        val layouts = layoutForest.buildReadLayouts(layoutReadCreateFunc)

        assertEquals(2, layouts.size)
        assertEquals("CAGCT", layouts.first().consensusSequenceString())
        assertEquals("11111", layouts.first().highQualSupportString())
        assertEquals("GCAGCT", layouts[1].consensusSequenceString())
        assertEquals("111111", layouts[1].highQualSupportString())
    }


    @Test
    fun testAddReadSimpleBranch()
    {
        // test very simple add high qual reads
        val layoutForest = LayoutForest(MIN_BASE_QUALITY, MIN_OVERLAP_BASES, MIN_SUPPORT_TO_SEAL_NODE)
        val seq1 = "CAGGT"
        val baseQual1 = SAMUtils.fastqToPhred("FFFFF") // F is 37, : is 25

        val read1 = LayoutForest.Read("read1", seq1.toByteArray(), baseQual1, 0)

        // seq 2 is different for the middle letter
        val seq2 = "CATGT"
        val baseQual2 = SAMUtils.fastqToPhred("FFFFF") // F is 37, : is 25
        val read2 = LayoutForest.Read("read2", seq2.toByteArray(), baseQual2, 0)

        // add it to the layout tree
        assertTrue(layoutForest.tryAddRead(read1))
        assertTrue(layoutForest.tryAddRead(read2))
        val layouts = layoutForest.buildReadLayouts(layoutReadCreateFunc).sortedBy({ l -> l.consensusSequenceString() })
        assertEquals(2, layouts.size)
        assertEquals(seq1, layouts[0].consensusSequenceString())
        assertEquals("11111", layouts[0].highQualSupportString())
        assertEquals(seq2, layouts[1].consensusSequenceString())
        assertEquals("11111", layouts[1].highQualSupportString())
    }

    @Test
    fun testAddReadLowQualMismatch()
    {
        // if we encounter a low qual mismatch, and can put it into other branch it is accepted

        val layoutForest = LayoutForest(MIN_BASE_QUALITY, MIN_OVERLAP_BASES, MIN_SUPPORT_TO_SEAL_NODE)
        val seq1 = "CAGGT"
        val baseQual1 = SAMUtils.fastqToPhred("FFFFF") // F is 37, : is 25

        val read1 = LayoutForest.Read("read1", seq1.toByteArray(), baseQual1, 0)

        // seq 2 is different for the middle letter
        val seq2 = "CATAT"
        val baseQual2 = SAMUtils.fastqToPhred("FF::F") // F is 37, : is 25
        val read2 = LayoutForest.Read("read2", seq2.toByteArray(), baseQual2, 0)

        // add it to the layout tree
        assertTrue(layoutForest.tryAddRead(read1))
        assertTrue(layoutForest.tryAddRead(read2))

        // since both match the high quality one we only get one layout
        val layouts = layoutForest.buildReadLayouts(layoutReadCreateFunc).sortedBy({ l -> l.consensusSequenceString() })
        assertEquals(1, layouts.size)
        assertEquals(seq1, layouts[0].consensusSequenceString())
        assertEquals("22112", layouts[0].highQualSupportString())
    }

    @Test
    fun testAddReadLowQualMismatch1()
    {
        // we first add the low qual one and then high qual one

        val layoutForest = LayoutForest(MIN_BASE_QUALITY, MIN_OVERLAP_BASES, MIN_SUPPORT_TO_SEAL_NODE)
        val seq1 = "CATAT"
        val baseQual1 = SAMUtils.fastqToPhred("FF::F") // F is 37, : is 25

        val read1 = LayoutForest.Read("read1", seq1.toByteArray(), baseQual1, 0)

        // seq 2 is different for the middle letter
        val seq2 = "CAGCT"
        val baseQual2 = SAMUtils.fastqToPhred("FFFFF") // F is 37, : is 25
        val read2 = LayoutForest.Read("read2", seq2.toByteArray(), baseQual2, 0)

        // add it to the layout tree
        assertTrue(layoutForest.tryAddRead(read1))
        assertTrue(layoutForest.tryAddRead(read2))

        // since both match the high quality one we only get one layout
        val layouts = layoutForest.buildReadLayouts(layoutReadCreateFunc).sortedBy({ l -> l.consensusSequenceString() })
        assertEquals(1, layouts.size)
        assertEquals(seq2, layouts[0].consensusSequenceString())
        assertEquals("22112", layouts[0].highQualSupportString())
    }

    @Test
    fun testLowQualBranch()
    {
        // we encounter a low qual branch that needs to be confirmed

        val layoutForest = LayoutForest(MIN_BASE_QUALITY, MIN_OVERLAP_BASES, MIN_SUPPORT_TO_SEAL_NODE)

        val seq1 = "CAGCT"
        val baseQual1 = SAMUtils.fastqToPhred("FFFFF") // F is 37, : is 25
        val read1 = LayoutForest.Read("read1", seq1.toByteArray(), baseQual1, 0)

        // seq 2 has a low qual diff in 3rd and high qual diff at 5th
        val seq2 = "CATCA"
        val baseQual2 = SAMUtils.fastqToPhred("FF:FF") // F is 37, : is 25
        val read2 = LayoutForest.Read("read2", seq2.toByteArray(), baseQual2, 0)

        // seq 3 is tricky, it actually flips the 3rd base to G, so the branch needs to be moved
        val seq3 = "CAGCA"
        val baseQual3 = SAMUtils.fastqToPhred("FFFFF") // F is 37, : is 25
        val read3 = LayoutForest.Read("read3", seq3.toByteArray(), baseQual3, 0)

        // add to the layout tree
        assertTrue(layoutForest.tryAddRead(read1))
        assertTrue(layoutForest.tryAddRead(read2))
        assertTrue(layoutForest.tryAddRead(read3))

        // since both match the high quality one we only get one layout
        val layouts = layoutForest.buildReadLayouts(layoutReadCreateFunc).sortedBy({ l -> l.consensusSequenceString() })
        assertEquals(2, layouts.size)
        assertEquals(seq3, layouts[0].consensusSequenceString())
        assertEquals("22122", layouts[0].highQualSupportString())
        assertEquals(seq1, layouts[1].consensusSequenceString())
        assertEquals("11111", layouts[1].highQualSupportString())
    }

    @Test
    fun testLowQualBranchExtend()
    {
        // we encounter a low qual branch that needs to be confirmed

        val layoutForest = LayoutForest(MIN_BASE_QUALITY, MIN_OVERLAP_BASES, MIN_SUPPORT_TO_SEAL_NODE)

        val seq1 = "CAGCT"
        val baseQual1 = SAMUtils.fastqToPhred("FF:FF") // F is 37, : is 25
        val read1 = LayoutForest.Read("read1", seq1.toByteArray(), baseQual1, 0)

        // seq 2 has a low qual diff in 3rd and high qual diff at 5th
        val seq2 = "CATCA"
        val baseQual2 = SAMUtils.fastqToPhred("FF:FF") // F is 37, : is 25
        val read2 = LayoutForest.Read("read2", seq2.toByteArray(), baseQual2, 0)

        // this one actually can merge with seq2, and extend it such that end result is CAGCAT
        val seq3 = "CAGTAT"
        val baseQual3 = SAMUtils.fastqToPhred("FFF:FF") // F is 37, : is 25
        val read3 = LayoutForest.Read("read3", seq3.toByteArray(), baseQual3, 0)

        // add to the layout tree
        assertTrue(layoutForest.tryAddRead(read1))
        assertTrue(layoutForest.tryAddRead(read2))
        assertTrue(layoutForest.tryAddRead(read3))

        // since both match the high quality one we only get one layout
        val layouts = layoutForest.buildReadLayouts(layoutReadCreateFunc).sortedBy({ l -> l.consensusSequenceString() })
        assertEquals(2, layouts.size)
        assertEquals("CAGCAT", layouts[0].consensusSequenceString())
        assertEquals("221121", layouts[0].highQualSupportString())
        assertEquals(seq1, layouts[1].consensusSequenceString())
        assertEquals("11011", layouts[1].highQualSupportString())
    }

    @Test
    fun testAddReadAfterError()
    {
        // test that layout tree allows reads to preferentially be added to most supported branch

        // we have following reads
        // GCAGCT
        //  CAGCT
        //   AGCTA  <--- last base is an error
        //   AGCTG
        //     CTGAA
        //      |
        //   aligned pos

        // resulting in following tree:
        // root
        //   \
        //    G1-C2-A4-G5-C5-T5-A1
        //                     \
        //                      G2-A1-A1
        // We want to short that all the reads except AGCTA should go the same layout.
        // There should be 2 layouts at the end, AGCTA and GCAGCTGAA

        val layoutForest = LayoutForest(MIN_BASE_QUALITY, MIN_OVERLAP_BASES, MIN_SUPPORT_TO_SEAL_NODE)

        layoutForest.tryAddRead(LayoutForest.Read("r1", 
            "GCAGCT".toByteArray(), SAMUtils.fastqToPhred("FFFFFF"), 5))

        layoutForest.tryAddRead(LayoutForest.Read("r2", 
            "CAGCT".toByteArray(), SAMUtils.fastqToPhred("FFFFF"), 4))

        layoutForest.tryAddRead(LayoutForest.Read("r3", 
            "AGCTA".toByteArray(), SAMUtils.fastqToPhred("FFFFF"), 3))

        layoutForest.tryAddRead(LayoutForest.Read("r4", 
            "AGCTG".toByteArray(), SAMUtils.fastqToPhred("FFFFF"), 3))

        layoutForest.tryAddRead(LayoutForest.Read("r5", 
            "CTGAA".toByteArray(), SAMUtils.fastqToPhred("FFFFF"), 1))

        // should have 2 sequences
        val layouts = layoutForest.buildReadLayouts(layoutReadCreateFunc).sortedByDescending({ l -> l.consensusSequenceString().length })
        assertEquals(2, layouts.size)
        assertEquals("GCAGCTGAA", layouts[0].consensusSequenceString())
        assertEquals("123344211", layouts[0].highQualSupportString())
        // 4 reads in first layout
        assertEquals(4, layouts[0].reads.size)
        assertEquals(5, layouts[0].alignedPosition)
        assertEquals("AGCTA", layouts[1].consensusSequenceString())
        assertEquals("11111", layouts[1].highQualSupportString())

        // 1 read in second layout
        assertEquals(1, layouts[1].reads.size)
        assertEquals(3, layouts[1].alignedPosition)
    }

    @Test
    fun testAddNewRootBranch()
    {
        // test that we can add a new branch to root correctly, if a read does not overlap existing branch

        // root
        //  \ \
        //   \ C1-A1-G1-C1-T1-G1-G1-A1-C1
        //    -------T1-C1-C1-A1

        val layoutForest = LayoutForest(MIN_BASE_QUALITY, MIN_OVERLAP_BASES, MIN_SUPPORT_TO_SEAL_NODE)

        val seq1 = "CAGCTGGAC"
        val baseQual1 = SAMUtils.fastqToPhred("F".repeat(seq1.length)) // F is 37, : is 25
        val read1 = LayoutForest.Read("read1", seq1.toByteArray(), baseQual1, 6)

        // read2 starts 2 bases after read1
        val seq2 = "TCCA"
        val baseQual2 = SAMUtils.fastqToPhred("FF:F") // F is 37, : is 25
        val read2 = LayoutForest.Read("read2", seq2.toByteArray(), baseQual2, 4)

        // add to the layout tree
        assertTrue(layoutForest.tryAddRead(read1))
        assertTrue(layoutForest.tryAddRead(read2))

        // should have 2 sequences
        var layouts = layoutForest.buildReadLayouts(layoutReadCreateFunc).sortedByDescending({ l -> l.consensusSequenceString().length })
        assertEquals(2, layouts.size)
        assertEquals(seq1, layouts[0].consensusSequenceString())
        assertEquals(6, layouts[0].alignedPosition)
        assertEquals(seq2, layouts[1].consensusSequenceString())
        assertEquals(4, layouts[1].alignedPosition)

        // test that we can extend the new branch
        val seq3 = "CCATT"
        val baseQual3 = SAMUtils.fastqToPhred("FFFFF") // F is 37, : is 25
        val read3 = LayoutForest.Read("read3", seq3.toByteArray(), baseQual3, 3)
        assertTrue(layoutForest.tryAddRead(read3))

        layouts = layoutForest.buildReadLayouts(layoutReadCreateFunc).sortedByDescending({ l -> l.consensusSequenceString().length })
        assertEquals(2, layouts.size)
        assertEquals(seq1, layouts[0].consensusSequenceString())
        assertEquals(6, layouts[0].alignedPosition)
        assertEquals("TCCATT", layouts[1].consensusSequenceString())
        assertEquals(4, layouts[1].alignedPosition)
    }

    @Test
    fun testAddReadToHighestSupport()
    {
        // test that we prefer to add read to branch with highest support, if both matches

        // A7-C7-T7-G7-G2-C1-C1-A1
        //            \
        //             T3-C3-C1-A1

        val layoutForest = LayoutForest(MIN_BASE_QUALITY, MIN_OVERLAP_BASES, MIN_SUPPORT_TO_SEAL_NODE)

        val seq1 = "CAGCT"
        val baseQual1 = SAMUtils.fastqToPhred("FF:FF") // F is 37, : is 25
        val read1 = LayoutForest.Read("read1", seq1.toByteArray(), baseQual1, 0)

        // seq 2 has a low qual diff in 3rd and high qual diff at 5th
        val seq2 = "CATCA"
        val baseQual2 = SAMUtils.fastqToPhred("FF:FF") // F is 37, : is 25
        val read2 = LayoutForest.Read("read2", seq2.toByteArray(), baseQual2, 0)

        // this one actually can merge with seq2, and extend it such that end result is CAGCAT
        val seq3 = "CAGTAT"
        val baseQual3 = SAMUtils.fastqToPhred("FFF:FF") // F is 37, : is 25
        val read3 = LayoutForest.Read("read3", seq3.toByteArray(), baseQual3, 0)

        // add to the layout tree
        assertTrue(layoutForest.tryAddRead(read1))
        assertTrue(layoutForest.tryAddRead(read2))
        assertTrue(layoutForest.tryAddRead(read3))

        // we get
        val layouts = layoutForest.buildReadLayouts(layoutReadCreateFunc).sortedBy({ l -> l.consensusSequenceString() })
        assertEquals(2, layouts.size)
        assertEquals("CAGCAT", layouts[0].consensusSequenceString())
        assertEquals("221121", layouts[0].highQualSupportString())
        assertEquals(seq1, layouts[1].consensusSequenceString())
        assertEquals("11011", layouts[1].highQualSupportString())
    }

    @Test
    fun testCannotAddToSealedNode()
    {
        // we have following reads
        // GCAGCT
        //  CAGCT
        //   AGCTGA
        //   AGCTGA
        //   AGCTAA  <---- this read cannot be added to the same branch since G is sealed
        //      |
        //   aligned pos

        // resulting in following tree:
        // root
        //   \
        //    G1-C2-A4-G4-C4-T4-G2-A2
        //   \
        //          A1-G1-C1-T1-A1-A1

        val layoutForest = LayoutForest(MIN_BASE_QUALITY, MIN_OVERLAP_BASES, 2)

        layoutForest.tryAddRead(LayoutForest.Read("r1",
            "GCAGCT".toByteArray(), SAMUtils.fastqToPhred("FFFFFF"), 5))

        layoutForest.tryAddRead(LayoutForest.Read("r2",
            "CAGCT".toByteArray(), SAMUtils.fastqToPhred("FFFFF"), 4))

        layoutForest.tryAddRead(LayoutForest.Read("r3",
            "AGCTGA".toByteArray(), SAMUtils.fastqToPhred("FFFFFF"), 3))

        layoutForest.tryAddRead(LayoutForest.Read("r4",
            "AGCTGA".toByteArray(), SAMUtils.fastqToPhred("FFFFFF"), 3))

        layoutForest.tryAddRead(LayoutForest.Read("r5",
            "AGCTAA".toByteArray(), SAMUtils.fastqToPhred("FFFFFF"), 3))

        // now want to check that the layoutForest root has two children
        assertEquals(2, layoutForest.roots.size)

        // should have 2 sequences
        val layouts = layoutForest.buildReadLayouts(layoutReadCreateFunc).sortedByDescending({ l -> l.consensusSequenceString().length })
        assertEquals(2, layouts.size)
        assertEquals("GCAGCTGA", layouts[0].consensusSequenceString())
        assertEquals("12444422", layouts[0].highQualSupportString())

        // 4 reads in first layout
        assertEquals(4, layouts[0].reads.size)
        assertEquals(5, layouts[0].alignedPosition)


        assertEquals("AGCTAA", layouts[1].consensusSequenceString())
        assertEquals("111111", layouts[1].highQualSupportString())

        // 1 read in second layout
        assertEquals(1, layouts[1].reads.size)
        assertEquals(3, layouts[1].alignedPosition)
    }

    @Test
    fun testReassignedReads()
    {
        // test that reads can be reassigned correctly
        // we have following reads
        // 1. GCAGCT
        // 2.   AGCTGAA   <---- this read needs to be reassigned
        // 3.   AGCTAA
        // 4.    GCTAAC
        // 5.  TAGCTGA
        //        |
        //     aligned pos

        // first 4 reads will go into one branch, read 5 will go to another
        // resulting in following tree:
        // root
        //   \
        //    G1-C1-A3-G4-C4-T4-G1-A1-A1
        //                     \A2-A2-C1
        //   \
        //    ---T1-A1-G1-C1-T1-G1-A1
        //
        // After first pass, read 3 will be reassigned to the 2nd branch
        // resulting in the following tree:
        // root
        //   \
        //    G1-C2-A4-G4-C4-T4-A2-A2-C1
        //   \
        //    ---T1-A2-G2-C2-T2-G2-A2-A1
        // we want to test that the reassignment occurred correctly

        val layoutForest = LayoutForest(MIN_BASE_QUALITY, MIN_OVERLAP_BASES, MIN_SUPPORT_TO_SEAL_NODE)

        layoutForest.tryAddRead(LayoutForest.Read("r1",
            "GCAGCT".toByteArray(), SAMUtils.fastqToPhred("FFFFFF"), 4))

        layoutForest.tryAddRead(LayoutForest.Read("r2",
            "AGCTGAA".toByteArray(), SAMUtils.fastqToPhred("FFFFFFF"), 2))

        layoutForest.tryAddRead(LayoutForest.Read("r3",
            "AGCTAA".toByteArray(), SAMUtils.fastqToPhred("FFFFFF"), 2))

        layoutForest.tryAddRead(LayoutForest.Read("r4",
            "GCTAAC".toByteArray(), SAMUtils.fastqToPhred("FFFFFF"), 1))

        layoutForest.tryAddRead(LayoutForest.Read("r5",
            "TAGCTGA".toByteArray(), SAMUtils.fastqToPhred("FFFFFFF"), 3))

        // now want to check that the layoutForest root has two children
        assertEquals(2, layoutForest.roots.size)

        // should have 2 sequences
        val layouts = layoutForest.buildReadLayouts(layoutReadCreateFunc).sortedByDescending({ l -> l.consensusSequenceString().length })
        assertEquals(2, layouts.size)
        assertEquals("GCAGCTAAC", layouts[0].consensusSequenceString())
        //assertEquals("12444422", layouts[0].highQualSupportString())

        // 3 reads in first layout
        assertEquals(3, layouts[0].reads.size)
        assertEquals(4, layouts[0].alignedPosition)


        assertEquals("TAGCTGAA", layouts[1].consensusSequenceString())
        assertEquals("12222221", layouts[1].highQualSupportString())

        // 2 read in second layout, after reassignment
        assertEquals(2, layouts[1].reads.size)
        assertEquals(3, layouts[1].alignedPosition)
    }

    companion object
    {
        const val MIN_BASE_QUALITY = 30.toByte()
        const val MIN_OVERLAP_BASES = 3
        const val MIN_SUPPORT_TO_SEAL_NODE = 3
    }
}
