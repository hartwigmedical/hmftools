package com.hartwig.hmftools.cider

import com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE
import com.hartwig.hmftools.common.bam.SamRecordUtils.getFivePrimeUnclippedPosition
import com.hartwig.hmftools.common.bam.SamRecordUtils.getThreePrimeUnclippedPosition
import com.hartwig.hmftools.common.bam.SamRecordUtils.mateUnmapped
import htsjdk.samtools.Cigar
import htsjdk.samtools.CigarElement
import htsjdk.samtools.CigarOperator
import htsjdk.samtools.SAMRecord
import org.apache.logging.log4j.LogManager
import java.io.InputStream
import java.util.*
import kotlin.io.path.createTempFile
import kotlin.math.abs
import kotlin.math.max

data class AlignmentBlock(val readStart: Int, val length: Int)

object CiderUtils
{
    private val sLogger = LogManager.getLogger(CiderUtils::class.java)

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

    // Calculate how much to trim a read to remove any adapter sequence leftover from sequencing.
    // This can occur if the fragment is very small, particularly for some types of panels.
    fun calculateAdapterSequenceTrim(read: SAMRecord): Pair<Int, Int>
    {
        // A = adapter DNA
        // C = align clip
        // fragment:             AAA--------------------------------AAA
        // left read (+):    5'     CCCC------------------------CCCCAAA  3'
        // right read (-):   3'  AAACCCC------------------------CCCC     5'
        // insert length:           --------------------------------
        // left trim:                                               ---
        // right trim:           ---

        if (!read.readUnmappedFlag && !mateUnmapped(read) && read.contig == read.mateReferenceName
            && read.readNegativeStrandFlag != read.mateNegativeStrandFlag)
        {
            val mateCigar = read.getStringAttribute(MATE_CIGAR_ATTRIBUTE)
            if (mateCigar != null)
            {
                val fivePrimePosition =
                    getFivePrimeUnclippedPosition(read.alignmentStart, read.cigarString, read.readNegativeStrandFlag)
                val threePrimePosition = getThreePrimeUnclippedPosition(read)
                val mateFivePrimePosition = getFivePrimeUnclippedPosition(
                    read.mateAlignmentStart, mateCigar, !read.mateNegativeStrandFlag
                )
                val insertLength = abs(mateFivePrimePosition - fivePrimePosition)
                val adapterBases = if (read.readNegativeStrandFlag)
                    max(0, mateFivePrimePosition - threePrimePosition) else
                    max(0, threePrimePosition - mateFivePrimePosition)
                if (insertLength < read.readLength && adapterBases > 0)
                {
                    return if (read.readNegativeStrandFlag) Pair(adapterBases, 0) else Pair(0, adapterBases)
                }
            }
        }

        // No adapter sequence or can't calculate it.
        return Pair(0, 0)
    }

    private const val ADAPTER_DNA_TRIM_ATTRIBUTE = "CiderDnaAdapterTrim"

    fun setAdapterDnaTrim(read: SAMRecord)
    {
        val trim = calculateAdapterSequenceTrim(read)
        // A bit ugly to store this into the SAMRecord, but it's better than replacing SAMRecord everywhere it's used.
        read.setTransientAttribute(ADAPTER_DNA_TRIM_ATTRIBUTE, trim)
        if (trim.first > 0 || trim.second > 0)
        {
            sLogger.trace("Trimmed adapter DNA {} {}:{}-{} trim={}",
                read.readName, read.referenceName, read.alignmentStart, read.alignmentEnd, trim)
        }
    }

    fun getAdapterDnaTrim(read: SAMRecord): Pair<Int, Int>
    {
        return (read.getTransientAttribute(ADAPTER_DNA_TRIM_ATTRIBUTE) ?: Pair(0, 0)) as Pair<Int, Int>
    }

    fun getResourceAsStream(name: String): InputStream
        = CiderUtils::class.java.classLoader.getResourceAsStream(name)!!

    fun getResourceAsFile(name: String, suffix: String? = null): String
    {
        // I don't think we can get a file directly from the resources, so have to write it to a temporary file.
        val resourceStream = getResourceAsStream(name)
        val file = createTempFile(suffix = suffix).toFile()
        file.deleteOnExit()
        file.writeBytes(resourceStream.readAllBytes())
        return file.path
    }
}