package com.hartwig.hmftools.cider

import com.hartwig.hmftools.common.codon.Nucleotides
import com.hartwig.hmftools.common.bam.SamRecordUtils
import com.hartwig.hmftools.common.utils.IntPair
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.util.SequenceUtil

// slice view of a read, we want to use this to help
// us manage times when we need to trim some poly G etc
// we do reverse complement first before slice if
// it is used
class ReadSlice(
    private val read: SAMRecord,
    val reverseComplement: Boolean,
    val sliceStart: Int,
    val sliceEnd: Int)
{
    init
    {
        require(sliceEnd > sliceStart, { "slice End <= sliceStart" })
    }

    val readName: String get()
    {
        return read.readName
    }

    val firstOfPairFlag: Boolean get()
    {
        return SamRecordUtils.firstInPair(read)
    }

    // just the data we want from here
    val readLength: Int get()
    {
        return sliceEnd - sliceStart
    }

    val readString: String get()
    {
        var readStr = read.readString
        if (reverseComplement)
            readStr = SequenceUtil.reverseComplement(readStr)
        return readStr.substring(sliceStart, sliceEnd)
    }

    val baseQualities: ByteArray get()
    {
        var bq = read.baseQualities
        if (reverseComplement)
            bq = bq.reversedArray()
        return bq.sliceArray(sliceStart until sliceEnd)
    }

    val baseQualityString: String get()
    {
        var bq = read.baseQualityString
        if (reverseComplement)
            bq = bq.reversed()
        return bq.substring(sliceStart, sliceEnd)
    }

    fun baseAt(slicePos: Int): Char
    {
        val b = read.readString[slicePositionToReadPosition(slicePos)]
        return if (reverseComplement) Nucleotides.complement(b) else b
    }

    fun baseQualityAt(slicePos: Int): Byte
    {
        return read.baseQualities[slicePositionToReadPosition(slicePos)]
    }

    fun readPositionToSlicePosition(readPos: Int): Int
    {
        val pos = if (reverseComplement) read.readLength - readPos - 1 else readPos
        return pos - sliceStart
    }

    fun slicePositionToReadPosition(slicePos: Int): Int
    {
        val pos = slicePos + sliceStart
        return if (reverseComplement) read.readLength - pos - 1 else pos
    }

    fun readRangeToSliceRange(readRangeStart: Int, readRangeEndExclusive: Int): IntPair
    {
        val start = if (reverseComplement) read.readLength - readRangeEndExclusive else readRangeStart
        val end = if (reverseComplement) read.readLength - readRangeStart else readRangeEndExclusive
        return IntPair(start - sliceStart, end - sliceStart)
    }

    fun sliceRangeToReadRange(sliceRangeStart: Int, sliceRangeEndExclusive: Int): IntPair
    {
        var start = sliceRangeStart + sliceStart
        var end = sliceRangeEndExclusive + sliceStart
        if (reverseComplement)
        {
            start = read.readLength - end
            end = read.readLength - start
        }
        return IntPair(start, end)
    }
}
