package com.hartwig.hmftools.cider

import htsjdk.samtools.Cigar
import htsjdk.samtools.CigarElement
import htsjdk.samtools.CigarOperator
import java.util.*

data class AlignmentBlock(val readStart: Int, val length: Int)

object CiderUtils
{
    fun calcBaseHash(base: Char): Int
    {
        return when (base)
        {
            'A' -> 0
            'T' -> 1
            'C' -> 2
            'G' -> 3
            else -> throw IllegalArgumentException("unknown base: $base")
        }
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
            if (pos < str.length && pos >= 0)
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
}