package com.hartwig.hmftools.cider

import htsjdk.samtools.Cigar
import htsjdk.samtools.CigarElement
import htsjdk.samtools.CigarOperator
import htsjdk.samtools.SAMRecord
import org.apache.logging.log4j.LogManager
import java.util.*

data class AlignmentBlock(val readStart: Int, val length: Int)

object CiderUtils
{
    private val sLogger = LogManager.getLogger(javaClass)

    fun conservedAA(vjGeneType: VJGeneType): Char
    {
        if (vjGeneType.vj == VJ.V)
            return 'C'
        if (vjGeneType == VJGeneType.IGHJ)
            return 'W'
        return 'F'
    }

    @JvmStatic
    fun countsToString(counts: IntArray): String
    {
        val RADIX = 36
        // create a support string
        val stringBuilder = StringBuilder(counts.size)

        for (count in counts)
        {
            if (count >= RADIX)
                stringBuilder.append('#')
            else
                stringBuilder.append(count.toString(radix = RADIX))
        }
        return stringBuilder.toString()
    }

    // insert dashes. Note that this works even if some dash positions
    // are invalid
    fun insertDashes(str: String, vararg dashPos: Int): String
    {
        val stringBuilder = StringBuilder(str)

        // we must go in reverse order
        for (i in dashPos.size - 1 downTo 0)
        {
            val pos = dashPos[i]
            if (pos < str.length && pos != 0)
            {
                stringBuilder.insert(pos, '-')
            }
        }

        return stringBuilder.toString()
    }

    // similar to SAMUtils.getAlignmentBlocks
    // but matched sections separated by I, D, P are merged together
    fun getAdjustedAlignmentBlocks(cigar: Cigar) : List<AlignmentBlock>
    {
        val alignmentBlocks: MutableList<AlignmentBlock> = ArrayList()
        var readBase = 1
        var currentBlockReadStart = -1
        for (e: CigarElement in cigar)
        {
            when (e.operator!!)
            {
                CigarOperator.P, CigarOperator.D ->
                {
                }
                CigarOperator.I -> readBase += e.length
                CigarOperator.M, CigarOperator.EQ, CigarOperator.X ->
                {
                    if (currentBlockReadStart == -1)
                    {
                        currentBlockReadStart = readBase
                    }
                    readBase += e.length
                }
                CigarOperator.S, CigarOperator.H, CigarOperator.N ->
                {
                    if (currentBlockReadStart != -1)
                    {
                        alignmentBlocks.add(AlignmentBlock(currentBlockReadStart, readBase - currentBlockReadStart))
                        currentBlockReadStart = -1
                    }
                    if (e.operator == CigarOperator.S)
                        readBase += e.length
                }
            }
        }
        if (currentBlockReadStart != -1)
        {
            alignmentBlocks.add(AlignmentBlock(currentBlockReadStart, readBase - currentBlockReadStart))
        }
        return Collections.unmodifiableList(alignmentBlocks)
    }

    //
    // apply reverse complement, trim bases and polyG trimming
    fun determineReadSlice(read: SAMRecord, useReverseComplement: Boolean, trimBases: Int) : ReadSlice?
    {
        // work out the slice start and end
        var sliceStart: Int = trimBases
        var sliceEnd: Int = read.readLength - trimBases

        // now we also want to try poly G tail trimming
        // we want to work out there the tail is.
        // the tail is on the right side and poly G if !read.readNegativeStrandFlag
        // the tail is on the left side and poly C otherwise
        if (!read.readNegativeStrandFlag)
        {
            // ends with poly G, but take trim bases into account
            val numGs = numTrailingPolyG(read.readString, sliceEnd)
            if (numGs >= CiderConstants.MIN_POLY_G_TRIM_COUNT)
            {
                sLogger.debug("read({}) strand(+) poly G tail of length({}) found({})",
                    read, numGs, read.readString)
                sliceEnd -= numGs + CiderConstants.POLY_G_TRIM_EXTRA_BASE_COUNT
            }
        }
        else
        {
            val numCs = numLeadingPolyC(read.readString, sliceStart)
            if (numCs >= CiderConstants.MIN_POLY_G_TRIM_COUNT)
            {
                sLogger.debug("read({}) strand(-) poly G tail of length({}) found({})",
                    read, numCs, read.readString)
                sliceStart += numCs + CiderConstants.POLY_G_TRIM_EXTRA_BASE_COUNT
            }
        }

        // the above logic is before reverse complement, the following logic is after
        // so we swap the start / end here
        if (useReverseComplement)
        {
            val sliceStartTmp = sliceStart
            sliceStart = read.readLength - sliceEnd
            sliceEnd = read.readLength - sliceStartTmp
        }

        if ((sliceEnd - sliceStart) < 5)
        {
            // if too little left don't bother
            return null
        }

        return ReadSlice(read, useReverseComplement, sliceStart, sliceEnd)
    }

    // apart from trim bases and poly G trim, this also remove bases that are "outside" the anchor
    // i.e. for V anchor read, bases before the anchor are trimmed, for J anchor read, bases after anchor are trimmed
    private fun determineAnchorReadSlice(read: SAMRecord, useReverseComplement: Boolean, vj: VJ,
                                   anchorOffsetStart: Int, anchorOffsetEnd: Int, trimBases: Int) : ReadSlice?
    {
        val readSlice = determineReadSlice(read, useReverseComplement, trimBases)

        if (readSlice == null)
            return null

        var sliceStart: Int = readSlice.sliceStart
        var sliceEnd: Int = readSlice.sliceEnd

        if (vj == VJ.V)
        {
            // for V, we layout from the anchor start, left to right
            // we are only interested in what comes after anchor start
            sliceStart = Math.max(anchorOffsetStart, sliceStart)
        }
        else
        {
            // for J, we layout from the anchor last, right to left
            // we are only interested in what comes before anchor end
            sliceEnd = Math.min(anchorOffsetEnd, sliceEnd)
        }

        if ((sliceEnd - sliceStart) < 5)
        {
            // if too little left don't bother
            return null
        }

        return ReadSlice(read, useReverseComplement, sliceStart, sliceEnd)
    }

    fun numTrailingPolyG(seq: String, sliceEnd: Int) : Int
    {
        for (i in 0 until sliceEnd)
        {
            if (seq[sliceEnd - i - 1] != 'G')
            {
                return i
            }
        }
        return sliceEnd
    }

    fun numLeadingPolyC(seq: String, sliceStart: Int) : Int
    {
        for (i in sliceStart until seq.length)
        {
            if (seq[i] != 'C')
            {
                return i - sliceStart
            }
        }
        return seq.length - sliceStart
    }

    fun safeSubstring(str: String, start: Int, end: Int): String
    {
        if (start >= str.length)
        {
            return ""
        }
        if (end >= str.length)
        {
            return str.substring(start)
        }
        return str.substring(start, end)
    }

    fun safeSubstring(str: String, range: IntRange): String
    {
        if (str.length <= range.first)
        {
            return ""
        }
        if (str.length <= range.last)
        {
            return str.substring(range.first)
        }
        return str.substring(range)
    }
}