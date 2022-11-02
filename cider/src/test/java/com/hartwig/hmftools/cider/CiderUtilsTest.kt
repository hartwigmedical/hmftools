package com.hartwig.hmftools.cider

import htsjdk.samtools.TextCigarCodec
import org.junit.Test
import kotlin.test.assertEquals

class CiderUtilsTest
{
    @Test
    fun testGetAdjustedAlignmentBlocks1()
    {
        var cigar = TextCigarCodec.decode("10S100M10S")!!
        var alignBlocks = CiderUtils.getAdjustedAlignmentBlocks(cigar)
        assertEquals(1, alignBlocks.size)
        assertEquals(AlignmentBlock(11, 100), alignBlocks[0])

        cigar = TextCigarCodec.decode("100M")!!
        alignBlocks = CiderUtils.getAdjustedAlignmentBlocks(cigar)
        assertEquals(1, alignBlocks.size)
        assertEquals(AlignmentBlock(1, 100), alignBlocks[0])

        cigar = TextCigarCodec.decode("100M10S")!!
        alignBlocks = CiderUtils.getAdjustedAlignmentBlocks(cigar)
        assertEquals(1, alignBlocks.size)
        assertEquals(AlignmentBlock(1, 100), alignBlocks[0])

        // small insert would merge together
        cigar = TextCigarCodec.decode("10S50M2I50M10S")!!
        alignBlocks = CiderUtils.getAdjustedAlignmentBlocks(cigar)
        assertEquals(1, alignBlocks.size)
        assertEquals(AlignmentBlock(11, 102), alignBlocks[0])

        // insert would merge together
        cigar = TextCigarCodec.decode("10S50M2I50M10S")!!
        alignBlocks = CiderUtils.getAdjustedAlignmentBlocks(cigar)
        assertEquals(1, alignBlocks.size)
        assertEquals(AlignmentBlock(11, 102), alignBlocks[0])

        // delete
        cigar = TextCigarCodec.decode("10S50M5D50M10S")!!
        alignBlocks = CiderUtils.getAdjustedAlignmentBlocks(cigar)
        assertEquals(1, alignBlocks.size)
        assertEquals(AlignmentBlock(11, 100), alignBlocks[0])

        // N, skipped section of ref genome
        cigar = TextCigarCodec.decode("10S50M1000N50M10S")!!
        alignBlocks = CiderUtils.getAdjustedAlignmentBlocks(cigar)
        assertEquals(2, alignBlocks.size)
        assertEquals(AlignmentBlock(11, 50), alignBlocks[0])
        assertEquals(AlignmentBlock(61, 50), alignBlocks[1])

        // hard clip
        cigar = TextCigarCodec.decode("10H50M1000N50M10H")!!
        alignBlocks = CiderUtils.getAdjustedAlignmentBlocks(cigar)
        assertEquals(2, alignBlocks.size)
        assertEquals(AlignmentBlock(1, 50), alignBlocks[0])
        assertEquals(AlignmentBlock(51, 50), alignBlocks[1])
    }
}