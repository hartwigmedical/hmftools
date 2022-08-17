package com.hartwig.hmftools.cider.layout

import com.hartwig.hmftools.cider.CiderUtils
import com.hartwig.hmftools.cider.ReadKey
import java.util.*
import kotlin.collections.ArrayList
import kotlin.collections.HashMap
import kotlin.collections.HashSet

//data class

class ReadLayout(var id: String = String())
{
    internal val allSequenceSupport: SequenceSupport = SequenceSupport()
    internal val highQualSequenceSupport: SequenceSupport = SequenceSupport()
    var alignedPosition: Int = 0

    data class Read (
        // have a source that allows us to refer back to where this comes from
        val source: Any,
        val readKey: ReadKey,
        val sequence: String,
        val baseQualities: ByteArray,
        val alignedPosition: Int)

    private val mutableReads: MutableList<Read> = ArrayList()
    private val mutableReadSet: MutableSet<ReadKey> = HashSet()
    val reads: List<Read> get() = mutableReads
    private var sequenceCache: String? = null

    val length: Int get() { return allSequenceSupport.support.size }
    val highQualSequence: String get() { return highQualSequenceSupport.sequence }

    fun consensusSequence(): String
    {
        if (sequenceCache == null)
            updateConsensusSequence()
        return sequenceCache!!
    }

    fun highQualSupportString(): String
    {
        return highQualSequenceSupport.supportString()
    }

    fun highQualSupportCounts(): IntArray
    {
        return highQualSequenceSupport.counts()
    }

    fun getAllSequenceSupportAt(index: Int) : Map.Entry<Char, Int>
    {
        return allSequenceSupport.support[index].likelyBaseSupport()
    }

    fun getHighQualSequenceSupportAt(index: Int) : Map.Entry<Char, Int>
    {
        return highQualSequenceSupport.support[index].likelyBaseSupport()
    }

    // get the reads covering the segment
    fun getSegmentReads(start: Int, end: Int) : List<Read>
    {
        val segmentReads = ArrayList<Read>()

        for (r in reads)
        {
            val readStart = getReadOffset(r)
            val readEnd = readStart + r.sequence.length

            if (readStart < end && start < readEnd)
            {
                segmentReads.add(r)
            }
        }
        return segmentReads
    }

    fun getReadsAt(index: Int) : List<Read>
    {
        return getSegmentReads(index, index + 1)
    }

    class BaseSupport
    {
        internal var likelyBase: Char = 'N'
        private var aCount: Int = 0
        private var cCount: Int = 0
        private var gCount: Int = 0
        private var tCount: Int = 0

        fun generateCountMap() : Map<Char, Int>
        {
            val countMap: MutableMap<Char, Int> = HashMap()

            for (base in arrayOf('A', 'C', 'G', 'T'))
            {
                val baseCount = count(base)

                if (baseCount > 0)
                    countMap[base] = baseCount
            }

            return countMap
        }

        fun count(base: Char): Int
        {
            return when (base)
                {
                    'A' -> aCount
                    'C' -> cCount
                    'G' -> gCount
                    'T' -> tCount
                    else -> 0
                }
        }

        // return true if the likely base has changed, false otherwise
        fun addToCount(base: Char) : Boolean
        {
            when (base)
            {
                'A' -> ++aCount
                'C' -> ++cCount
                'G' -> ++gCount
                'T' -> ++tCount
            }
            return updateLikelyBase(base)
        }

        private fun updateLikelyBase(baseWithAddedCount: Char) : Boolean
        {
            if (likelyBase != baseWithAddedCount &&
                count(baseWithAddedCount) > count(likelyBase))
            {
                likelyBase = baseWithAddedCount
                return true
            }
            return false
        }

        fun likelyBaseSupport(): Map.Entry<Char, Int>
        {
            if (likelyBase == 'N')
                return sNullSupport
            return AbstractMap.SimpleImmutableEntry(likelyBase, likelyBaseSupportCount())
        }

        fun likelyBaseSupportCount(): Int
        {
            return count(likelyBase)
        }
    }

    class SequenceSupport (
        var sequence: String = String())
    {
        private val mMutableSupport: MutableList<BaseSupport> = ArrayList()
        val support: List<BaseSupport> get() { return mMutableSupport }

        internal fun ensureLength(l: Int)
        {
            val lengthChange = l - support.size
            if (lengthChange > 0)
            {
                for (i in 0 until lengthChange)
                {
                    mMutableSupport.add(BaseSupport())
                }
                sequence += "N".repeat(lengthChange)
            }
        }

        internal fun shiftRightBy(s: Int)
        {
            if (s <= 0)
                return
            // also need to add some items
            for (i in 0 until s)
            {
                mMutableSupport.add(BaseSupport())
            }
            // rotate the new items to the front
            Collections.rotate(support, s)
            sequence = "N".repeat(s) + sequence
        }

        internal fun updateSequence()
        {
            val stringBuilder = StringBuilder(support.size)
            for (baseSupport in support)
            {
                stringBuilder.append(baseSupport.likelyBase)
            }
            sequence = stringBuilder.toString()
        }

        fun supportString(): String
        {
            return CiderUtils.countsToString(counts())
        }

        fun counts(): IntArray
        {
            val countArray = IntArray(sequence.length)

            for (i in sequence.indices)
            {
                val base: Char = sequence[i]
                val readCount: Int = support[i].count(base)
                countArray[i] = readCount
            }
            return countArray
        }
    }

    fun addRead(read: Read, minBaseQuality: Byte) : Boolean
    {
        if (!mutableReadSet.add(read.readKey))
            return false

        if (mutableReads.isEmpty())
        {
            // this if first read
            alignedPosition = read.alignedPosition
        }
        else if (read.alignedPosition > alignedPosition)
        {
            // update max offset for each group item
            val offsetChange = read.alignedPosition - alignedPosition

            if (offsetChange > 0)
            {
                allSequenceSupport.shiftRightBy(offsetChange)
                highQualSequenceSupport.shiftRightBy(offsetChange)
            }

            alignedPosition = read.alignedPosition
        }

        mutableReads.add(read)

        // now use the aligned position to line up the read with the group
        val seqOffset = getReadOffset(read)
        val readSeqEnd = read.sequence.length + seqOffset

        allSequenceSupport.ensureLength(readSeqEnd)
        highQualSequenceSupport.ensureLength(readSeqEnd)

        // now we want to update the consensus sequence
        // we recalculate it from scratch for now

        var allSequenceChanged = false
        var highQualSequenceChanged = false

        // update the support
        for (i in read.sequence.indices)
        {
            val baseQual = read.baseQualities[i]

            val b: Char = read.sequence[i]
            if (b != 'N')
            {
                allSequenceChanged = allSequenceSupport.support[seqOffset + i].addToCount(b) or allSequenceChanged

                if (baseQual >= minBaseQuality)
                {
                    highQualSequenceChanged = highQualSequenceSupport.support[seqOffset + i].addToCount(b) or highQualSequenceChanged
                }
            }
        }

        if (allSequenceChanged)
            allSequenceSupport.updateSequence()
        if (highQualSequenceChanged)
            highQualSequenceSupport.updateSequence()

        if (allSequenceChanged || highQualSequenceChanged)
            sequenceCache = null

        return true
    }

    private fun updateConsensusSequence()
    {
        // combine high qual with low qual
        val stringBuilder = StringBuilder(allSequenceSupport.sequence.length)
        for (i in allSequenceSupport.sequence.indices)
        {
            val base = highQualSequenceSupport.sequence[i]
            if (base == 'N')
            {
                // if high quality sequence is not sure, use the low qual one
                stringBuilder.append(allSequenceSupport.sequence[i])
            }
            else
            {
                stringBuilder.append(base)
            }
        }
        sequenceCache = stringBuilder.toString()
    }

    private fun getReadOffset(read: Read) : Int
    {
        return alignedPosition - read.alignedPosition
    }

    // merge all the reads from another layout. alignedPositionShift specifies the change to
    // the aligned position of the reads that are merged in.
    fun mergeIn(layout: ReadLayout, alignedPositionShift: Int, minBaseQuality: Byte)
    {
        for (r in layout.reads)
        {
            if (alignedPositionShift == 0)
                addRead(r, minBaseQuality)
            else
                addRead(r.copy(alignedPosition = r.alignedPosition + alignedPositionShift), minBaseQuality)
        }
    }

    companion object
    {
        val sNullSupport : Map.Entry<Char, Int> = AbstractMap.SimpleImmutableEntry('N', 0)

        // merge in another layout to this layout, and create a new layout
        // the new layout will have the same aligned position as this layout
        fun merge(layout1: ReadLayout, layout2: ReadLayout,
                  alignedPositionShift1: Int, alignedPositionShift2: Int,
                  minBaseQuality: Byte) : ReadLayout
        {
            // unfortunately in java / kotlin it is difficult to clone objects
            // will do it the slow way of building it from ground up again
            val merged = ReadLayout()

            merged.mergeIn(layout1, alignedPositionShift1, minBaseQuality)
            merged.mergeIn(layout2, alignedPositionShift2, minBaseQuality)
            return merged
        }
    }
}
