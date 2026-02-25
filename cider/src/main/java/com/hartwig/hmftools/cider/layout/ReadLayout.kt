package com.hartwig.hmftools.cider.layout

import com.hartwig.hmftools.cider.CiderUtils
import com.hartwig.hmftools.cider.ReadKey
import htsjdk.samtools.util.SequenceUtil.N
import htsjdk.samtools.util.SequenceUtil.A
import htsjdk.samtools.util.SequenceUtil.T
import htsjdk.samtools.util.SequenceUtil.C
import htsjdk.samtools.util.SequenceUtil.G
import java.util.*
import kotlin.collections.ArrayList
import kotlin.collections.HashSet

open class ReadLayout(var id: String = String())
{
    internal val allSequenceSupport: SequenceSupport = SequenceSupport()
    internal val highQualSequenceSupport: SequenceSupport = SequenceSupport()
    var alignedPosition: Int = 0

    abstract class Read (
        val readKey: ReadKey,
        val sequence: ByteArray,
        val baseQualities: ByteArray,
        val alignedPosition: Int)
    {
        init
        {
            require(sequence.size == baseQualities.size, { "sequence.length != baseQualities.size" })
        }

        val readLength: Int get() { return sequence.size }
        val sequenceString: String get() { return String(sequence) }

        // this function allows us to copy the layout read with different aligned position
        abstract fun copy(alignedPosition: Int): Read
    }

    private val mutableReads: MutableList<Read> = ArrayList()
    private val mutableReadSet: MutableSet<ReadKey> = HashSet()
    val reads: List<Read> get() = mutableReads
    private var sequenceCache: ByteArray? = null

    val length: Int get() { return allSequenceSupport.support.size }

    // high qual sequence would contain N if there is no high quality base at a position
    // val highQualSequence: String get() { return highQualSequenceSupport.sequence }

    // consensus sequence is high qual sequence but N replace with
    // the low quality sequence
    fun consensusSequence(): ByteArray
    {
        if (sequenceCache == null)
            updateConsensusSequence()
        return sequenceCache!!
    }

    fun consensusSequenceString(): String
    {
        return String(consensusSequence())
    }

    //
    // highQualReadFraction is the % number of high quality reads supports each base
    // i.e. if there are 1000 reads in this layout, and highQualReadFraction = 0.01
    //      then each base needs to be supported by 10 high quality reads or it will be set
    //      to N
    fun highConfidenceSequence(highQualReadFraction: Double): ByteArray
    {
        val minHighQualCount = Math.max((reads.size * highQualReadFraction).toInt(), 1)
        val highQualCounts = highQualSupportCounts()
        val seq = ByteArray(highQualSequenceSupport.sequence.size)
        for (i in highQualSequenceSupport.sequence.indices)
        {
            val base = highQualSequenceSupport.sequence[i]

            if (highQualCounts[i] >= minHighQualCount)
            {
                seq[i] = base
            }
            else
            {
                seq[i] = N
            }
        }
        return seq
    }

    fun highQualSupportString(): String
    {
        return highQualSequenceSupport.supportString()
    }

    fun highQualSupportCounts(): IntArray
    {
        return highQualSequenceSupport.counts()
    }

    fun getHighQualSequenceSupportAt(index: Int) : Map.Entry<Byte, Int>
    {
        return highQualSequenceSupport.support[index].likelyBaseSupport()
    }

    class BaseSupport
    {
        internal var likelyBase: Byte = N
        private var aCount: Int = 0
        private var cCount: Int = 0
        private var gCount: Int = 0
        private var tCount: Int = 0

        fun count(base: Byte): Int
        {
            return when (base)
                {
                    A -> aCount
                    C -> cCount
                    G -> gCount
                    T -> tCount
                    else -> 0
                }
        }

        // return true if the likely base has changed, false otherwise
        fun addToCount(base: Byte) : Boolean
        {
            when (base)
            {
                A -> ++aCount
                C -> ++cCount
                G -> ++gCount
                T -> ++tCount
            }
            return updateLikelyBase(base)
        }

        private fun updateLikelyBase(baseWithAddedCount: Byte) : Boolean
        {
            if (likelyBase != baseWithAddedCount &&
                count(baseWithAddedCount) > count(likelyBase))
            {
                likelyBase = baseWithAddedCount
                return true
            }
            return false
        }

        fun likelyBaseSupport(): Map.Entry<Byte, Int>
        {
            if (likelyBase == N)
                return sNullSupport
            return AbstractMap.SimpleImmutableEntry(likelyBase, likelyBaseSupportCount())
        }

        fun likelyBaseSupportCount(): Int
        {
            return count(likelyBase)
        }
    }

    class SequenceSupport (
        var sequence: ByteArray = ByteArray(0))
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
                sequence += ByteArray(lengthChange) { N }
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
            sequence = ByteArray(s) { N } + sequence
        }

        internal fun updateSequence()
        {
            if (sequence.size != support.size)
                sequence = ByteArray(support.size)

            for (i in support.indices)
            {
                sequence[i] = support[i].likelyBase
            }
        }

        fun supportString(): String
        {
            return CiderUtils.countsToString(counts())
        }

        fun counts(): IntArray
        {
            val countArray = IntArray(sequence.size)

            for (i in sequence.indices)
            {
                countArray[i] = support[i].count(sequence[i])
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
        val readSeqEnd = read.sequence.size + seqOffset

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

            val b: Byte = read.sequence[i]
            if (b != N)
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
        val size = length
        if (sequenceCache == null || sequenceCache!!.size != size)
            sequenceCache = ByteArray(size)

        // combine high qual with low qual
        for (i in allSequenceSupport.sequence.indices)
        {
            val base = highQualSequenceSupport.sequence[i]
            if (base == N)
            {
                // if high quality sequence is not sure, use the low qual one
                sequenceCache!![i] = allSequenceSupport.sequence[i]
            }
            else
            {
                sequenceCache!![i] = base
            }
        }
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
        val sNullSupport : Map.Entry<Byte, Int> = AbstractMap.SimpleImmutableEntry(N, 0)

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
