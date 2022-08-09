package com.hartwig.hmftools.cdr3.layout

import org.apache.logging.log4j.LogManager
import java.util.*
import kotlin.collections.ArrayList

class ReadLayoutBuilder(inputReads: List<ReadLayout.Read>, minBaseQuality: Int, minMatchedBases: Int, minMatchRatio: Double, alignLeft: Boolean)
{
    val minBaseQuality: Byte = minBaseQuality.toByte()
    val minMatchedBases: Int = minMatchedBases
    val minMatchRatio: Double = minMatchRatio
    val inputReadList: List<ReadLayout.Read>

    init
    {
        val readDataMutableList = ArrayList<ReadLayout.Read>()
        readDataMutableList.addAll(inputReads)

        if (alignLeft)
        {
            // we want to cluster them based on aligned position so the sequence is built up from
            // longest to shortest, highest quality first
            readDataMutableList.sortWith(Collections.reverseOrder(
                    Comparator.comparingInt({ r: ReadLayout.Read -> r.alignedPosition + r.sequence.length })
                        .thenComparingDouble({ r: ReadLayout.Read -> r.baseQualities.average() }) // handle the highest quality ones first
                        .thenComparingInt({ r: ReadLayout.Read -> r.alignedPosition })
                        .thenComparing({ r: ReadLayout.Read -> r.readKey.readName }) // lastly we use read Id just in case
            ))
        }
        else
        {
            // start with sequences that are longest, highest quality first
            readDataMutableList.sortWith(Collections.reverseOrder(
                Comparator.comparingInt({ r: ReadLayout.Read -> r.alignedPosition })
                    .thenComparingDouble({ r: ReadLayout.Read -> r.baseQualities.average() }) // handle the highest quality ones first
                    .thenComparingInt({ r: ReadLayout.Read -> r.alignedPosition + r.sequence.length })
                    .thenComparing({ r: ReadLayout.Read -> r.readKey.readName }) // lastly we use read Id just in case
            ))
        }

        inputReadList = readDataMutableList
    }

    fun build(): List<ReadLayout>
    {
        sLogger.info("building layouts from {} reads", inputReadList.size)

        val readLayouts: MutableList<ReadLayout> = ArrayList()

        // go through the read data list, and add one by one to the list of clusters
        // if there are multiple clusters that matches, we choose the highest one
        for (readData in inputReadList)
        {
            var matchedLayout: ReadLayout? = null

            for (i in readLayouts.indices)
            {
                val readLayout = readLayouts[i]

                if (matchedLayout != null && matchedLayout.reads.size > readLayout.reads.size)
                {
                    // no need to try this one
                    continue
                }

                // test the read against this overlay
                if (layoutMatch(readLayout, readData, minBaseQuality, minMatchedBases))
                {
                    //sLogger.info("match found")
                    matchedLayout = readLayout
                }
            }

            if (matchedLayout == null)
            {
                val layout = ReadLayout()
                layout.addRead(readData, minBaseQuality)
                readLayouts.add(layout)
            }
            else
            {
                // add to all layouts that match
                // for (l in matchedLayouts)
                    //addToOverlay(l, readData, minBaseQuality)
                matchedLayout.addRead(readData, minBaseQuality)
            }
        }

        sLogger.info("built {} layouts from {} reads", readLayouts.size, inputReadList.size)

        return readLayouts
    }

    companion object
    {
        private val sLogger = LogManager.getLogger(ReadLayoutBuilder::class.java)

        //
        @JvmStatic
        fun layoutMatch(layout: ReadLayout, readData: ReadLayout.Read, baseQualityCutoff: Byte, minOverlap: Int): Boolean
        {
            // NOTE: alignedPosition is with respect to itself.
            // i.e. for a layout it is the position from layout start
            // for a read it is position from the read sequence start
            // it can actually be negative or larger than the sequence

            val alignedPosMin = Math.min(layout.alignedPosition, readData.alignedPosition)

            // this will work also if the position is negative, i.e. before the sequence
            val layoutOffsetStart: Int = layout.alignedPosition - alignedPosMin
            val readOffsetStart: Int = readData.alignedPosition - alignedPosMin

            // n is the number of bases overlap
            val n: Int = Math.min(layout.length - layoutOffsetStart, readData.sequence.length - readOffsetStart)

            if (n < minOverlap)
                return false

            // we only compare the CDR3 DNA part, the other part we ignore
            return sequenceMatch(layout.highQualSequence, readData.sequence,
                null, readData.baseQualities,
                layoutOffsetStart, readOffsetStart,
                n,
                baseQualityCutoff)
        }

        // if base qual is null then don't check it
        @JvmStatic
        fun sequenceMatch(read1Seq: String, read2Seq: String,
                               read1Qual: ByteArray?, read2Qual: ByteArray?,
                               read1StartOffset: Int, read2StartOffset: Int,
                               n: Int,
                               baseQualityCutoff: Byte): Boolean
        {
            var i1: Int = read1StartOffset
            var i2: Int = read2StartOffset
            var comparedCount: Int = 0
            var matchCount: Int = 0

            if (read1StartOffset + n > read1Seq.length)
            {
                throw RuntimeException("read1StartOffset + n > read1Seq.length")
            }

            if (read2StartOffset + n > read2Seq.length)
            {
                throw RuntimeException("read2StartOffset + n > read2Seq.length")
            }

            for (i in 0 until n)
            {
                val b1: Char = read1Seq[i1]
                val b2: Char = read2Seq[i2]

                if (b1 == b2)
                {
                    // TODO: check with Peter
                    ++comparedCount
                    ++matchCount
                }
                else if ((read1Qual == null || read1Qual[i1] >= baseQualityCutoff) &&
                    (read2Qual == null || read2Qual[i2] >= baseQualityCutoff) &&
                    b1 != 'N' && b2 != 'N')
                {
                    // if both high quality and no N then we say this is a real mismatch
                    ++comparedCount
                    return false
                }

                ++i1
                ++i2
            }

            return true
        }

        // if base qual is null then don't check it
        @JvmStatic
        fun sequenceMatchCount(read1Seq: String, read2Seq: String,
                               read1Qual: ByteArray?, read2Qual: ByteArray?,
                               read1StartOffset: Int, read2StartOffset: Int,
                               n: Int,
                               baseQualityCutoff: Byte): Pair<Int, Int>
        {
            var i1: Int = read1StartOffset
            var i2: Int = read2StartOffset
            var comparedCount: Int = 0
            var matchCount: Int = 0

            if (read1StartOffset + n > read1Seq.length)
            {
                throw RuntimeException("read1StartOffset + n > read1Seq.length")
            }

            if (read2StartOffset + n > read2Seq.length)
            {
                throw RuntimeException("read2StartOffset + n > read2Seq.length")
            }

            for (i in 0 until n)
            {
                val b1: Char = read1Seq[i1]
                val b2: Char = read2Seq[i2]

                if (b1 == b2)
                {
                    // TODO: check with Peter
                    ++comparedCount
                    ++matchCount
                }
                else if ((read1Qual == null || read1Qual[i1] >= baseQualityCutoff) &&
                        (read2Qual == null || read2Qual[i2] >= baseQualityCutoff) &&
                        b1 != 'N' && b2 != 'N')
                {
                    // if both high quality and no N then we say this is a real mismatch
                    ++comparedCount
                }

                ++i1
                ++i2
            }

            return Pair(matchCount, comparedCount)
        }
    }
}
