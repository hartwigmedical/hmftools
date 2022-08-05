package com.hartwig.hmftools.cdr3.layout

import com.hartwig.hmftools.common.utils.Doubles
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
                //val (matchCount, compareCount) = layoutMatch(readLayout, readData, alignLeft, minBaseQuality)
                //val matchRatio = matchCount / compareCount.toDouble()

                if (layoutMatch(readLayout, readData, minBaseQuality, minMatchedBases))
                {
                    //sLogger.info("match found")
                    matchedLayout = readLayout
                }
            }

            if (matchedLayout == null)
            {
                val overlay = ReadLayout()
                addToOverlay(overlay, readData, minBaseQuality)
                readLayouts.add(overlay)
            }
            else
            {
                // add to all layouts that match
                // for (l in matchedLayouts)
                    //addToOverlay(l, readData, minBaseQuality)
                addToOverlay(matchedLayout, readData, minBaseQuality)
            }
        }

        sLogger.info("built {} layouts from {} reads", readLayouts.size, inputReadList.size)

        return readLayouts
    }

    /*
    fun clusterReads(readDataList: ArrayList<VJReadData>, mismatchCutoff: Double, minMatchedBases: Int, minBaseQuality: Byte)
    {
        // create a list of all pairs
        val matrix = Matrix(readDataList.size, readDataList.size)

        for (i in 0 until readDataList.size)
        {
            val read1 = readDataList[i]
            assert(read1.index == i)

            for (j in i + 1 until readDataList.size)
            {
                val read2 = readDataList[j]
                assert(read2.index == j)

                val similarity = scoreReadPair(read1, read2, minBaseQuality, minMatchedBases)

                matrix[i, j] = similarity
            }
        }

        // now we need to define a direction to process each of these pairs
        // we should process the longest reads
        // and the most similar reads
    }*/

    companion object
    {
        private val sLogger = LogManager.getLogger(ReadLayoutBuilder::class.java)

        fun addToOverlay(group: ReadLayout, read: ReadLayout.Read, minBaseQuality: Byte)
        {
            group.addRead(read, minBaseQuality)

            /*

            if (readData.alignedPosition > group.alignedPosition)
            {
                // update max offset for each group item
                val offsetChange = readData.alignedPosition - group.alignedPosition
                group.reads.forEach({ o -> o.overlayOffset += offsetChange })

                // also need to add some items
                for (i in 0 until offsetChange)
                {
                    group.support.add(BaseSupport())
                }

                // rotate the new items to the front
                Collections.rotate(group.support, offsetChange)

                group.alignedPosition = readData.alignedPosition
            }

            // now use the aligned position to line up the read with the group
            val seqOffset = group.alignedPosition - readData.alignedPosition

            group.reads.add(
                ReadLayout.Read(
                    readData.key, readData.sequence, readData.baseQualities,
                    seqOffset
                )
            )

            for (i in group.support.size until readData.sequence.length + seqOffset)
            {
                group.support.add(BaseSupport())
            }

            // now we want to update the consensus sequence
            // we recalculate it from scratch for now

            // update the support
            for (i in readData.sequence.indices)
            {
                val baseQual = readData.baseQualities[i]

                if (baseQual < minBaseQuality)
                    continue

                val b: Char = readData.sequence[i]
                if (b != 'N')
                    group.support[seqOffset + i].addToCount(b)
            }

            // now read out the sequence
            val stringBuilder = StringBuilder(group.support.size)

            for (baseSupport in group.support)
            {
                stringBuilder.append(baseSupport.likelyBase())
            }

            group.sequence = stringBuilder.toString()

             */
        }

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
