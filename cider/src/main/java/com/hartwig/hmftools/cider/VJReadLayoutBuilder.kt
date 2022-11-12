package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.layout.ReadLayout
import htsjdk.samtools.SAMRecord
import org.apache.logging.log4j.LogManager
import java.util.*
import kotlin.collections.ArrayList

// helper class to convert from the outer VJ classes to the layout classes
// create an interface to make it easier to test
interface IVJReadLayoutAdaptor
{
    //fun toReadCandidate(read: ReadLayout.Read) : VJReadCandidate
    fun getAnchorMatchMethod(layout: ReadLayout) : VJReadCandidate.MatchMethod
    fun getTemplateAnchorSequence(layout: ReadLayout) : String
    fun getAnchorRange(vj: VJ, layout: ReadLayout) : IntRange?
}

// This class is the mapper between the candidate reads and the layout object
// There are several hacks that is used here to keep track of
//
// 1. for layout that comes from V read candidates, the aligned position is the last base of the V anchor i.e.
//                                  * <-- this T is the aligned position of the layout
//    AGATCTGAG-GACACGGCCGTGTATTACTGT-GCGAGAGACACAGTGTGAAAACCCACATCCTGAGAGTGTCAGAAACCCTGAGGGA
//              |___________________|
//                V anchor
// 2. for layout that comes from J read candidates, the aligned position is the first base of the J anchor, i.e.
//         this C is the aligned position of the layout --> *
//    AGATCTGAGGACACGGCCGTGTATTACTGTGCGAGAGACACAGTGTGAAAACC-CACATCCTGAGAGTGTCAGAA-ACCCTGAGGGA
//                                                          |___________________|
//                                                                 J anchor
// Most functions here rely on this.
class VJReadLayoutBuilder(private val trimBases: Int, private val minBaseQuality: Int) : IVJReadLayoutAdaptor
{
    private class VjLayoutRead private constructor(
        val layoutReadSlice: ReadSlice,
        val readCandidate: VJReadCandidate,
        readKey: ReadKey,
        sequence: String,
        baseQualities: ByteArray,
        alignedPosition: Int)
        : ReadLayout.Read(readKey, sequence, baseQualities, alignedPosition)
    {
        constructor(layoutReadSlice: ReadSlice, readCandidate: VJReadCandidate, alignedPosition: Int)
                : this(layoutReadSlice,
            readCandidate,
            ReadKey(layoutReadSlice.readName, layoutReadSlice.firstOfPairFlag),
            layoutReadSlice.readString,
            layoutReadSlice.baseQualities,
            alignedPosition)
        {
        }

        override fun copy(alignedPosition: Int): ReadLayout.Read
        {
            return VjLayoutRead(layoutReadSlice, readCandidate, readKey, sequence, baseQualities, alignedPosition)
        }
    }

    fun readCandidateToLayoutRead(readCandidate: VJReadCandidate) : ReadLayout.Read?
    {
        val slice = determineReadSlice(readCandidate.read, readCandidate.useReverseComplement)

        if (slice == null)
            return null

        // now determine the aligned position
        val alignedPosition: Int

        if (readCandidate.vjGeneType.vj == VJ.V)
        {
            // aligned position we must take into account that we remove all bases before vAnchor start
            alignedPosition = readCandidate.anchorOffsetEnd - 1 - slice.sliceStart
        }
        else
        {
            // aligned position we must take into account that we remove all bases before vAnchor start
            alignedPosition = readCandidate.anchorOffsetStart - slice.sliceStart
        }

        return VjLayoutRead(slice, readCandidate, alignedPosition)
    }

    fun toReadCandidate(read: ReadLayout.Read) : VJReadCandidate
    {
        return (read as VjLayoutRead).readCandidate
    }

    fun toLayoutReadSlice(read: ReadLayout.Read) : ReadSlice
    {
        return (read as VjLayoutRead).layoutReadSlice
    }

    fun getReadCandidates(layout: ReadLayout) : List<VJReadCandidate>
    {
        return layout.reads.map({ read: ReadLayout.Read -> toReadCandidate(read)})
    }

    // get the anchor boundary position for this layout read
    fun getAnchorBoundaryPosition(read: ReadLayout.Read) : Int
    {
        return read.alignedPosition
    }

    // apply trim bases and polyG trimming
    private fun determineReadSlice(read: SAMRecord, useReverseComplement: Boolean) : ReadSlice?
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
            val numGs = CiderUtils.numTrailingPolyG(read.readString, sliceEnd)
            if (numGs >= CiderConstants.MIN_POLY_G_TRIM_COUNT)
            {
                sLogger.debug("read({}) strand(+) poly G tail of length({}) found({})",
                    read, numGs, read.readString)
                sliceEnd -= numGs + CiderConstants.POLY_G_TRIM_EXTRA_BASE_COUNT
            }
        }
        else
        {
            val numCs = CiderUtils.numLeadingPolyC(read.readString, sliceStart)
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

    override fun getAnchorMatchMethod(layout: ReadLayout): VJReadCandidate.MatchMethod
    {
        val readCandidates = getReadCandidates(layout)

        // we need to get a few values from read candidates
        if (readCandidates.isEmpty())
            throw IllegalArgumentException("read candidate list is empty")

        // just return first one for now, not the best but should be fine
        return readCandidates.first().matchMethod
    }

    override fun getTemplateAnchorSequence(layout: ReadLayout) : String
    {
        val readCandidates = getReadCandidates(layout)

        // we need to get a few values from read candidates
        if (readCandidates.isEmpty())
            throw IllegalArgumentException("read candidate list is empty")

        return readCandidates.maxByOrNull({ r -> r.similarityScore })?.templateAnchorSequence ?: ""
    }

    override fun getAnchorRange(vj: VJ, layout: ReadLayout) : IntRange?
    {
        val layoutReads = layout.reads.map { o: ReadLayout.Read -> toReadCandidate(o) }
            .toList()

        // we more or less get the top one
        val anchorLength = layoutReads.maxOfOrNull { o: VJReadCandidate -> o.anchorOffsetEnd - o.anchorOffsetStart } ?: 0

        val anchorRange =
        // for V read we align to last base of anchor, for J read we align to first base of the anchor
        if (vj == VJ.V)
            layout.alignedPosition - anchorLength + 1 .. layout.alignedPosition
        else if (vj == VJ.J)
            layout.alignedPosition until layout.alignedPosition + anchorLength
        else
            return null

        if (anchorRange.first >= layout.length || anchorRange.last < 0)
            return null

        // protect against 0 and end
        return Math.max(0, anchorRange.first) until Math.min(anchorRange.last + 1, layout.length)
    }

    fun getAnchorRange(geneType: VJGeneType, layout: ReadLayout) : IntRange?
    {
        return getAnchorRange(geneType.vj, layout)
    }

    fun buildLayouts(geneType: VJGeneType, readCandidates: List<VJReadCandidate>,
                     minMatchedBases: Int, minMatchRatio: Double)
    : List<ReadLayout>
    {
        sLogger.info("building {} layouts from {} reads", geneType, readCandidates.size)

        val layoutReads = ArrayList<ReadLayout.Read>()

        for (r in readCandidates)
        {
            val layoutRead = readCandidateToLayoutRead(r)
            if (layoutRead != null)
                layoutReads.add(layoutRead)
        }

        if (geneType.vj == VJ.V)
        {
            // we want to cluster them based on aligned position so the sequence is built up from
            // longest to shortest, highest quality first
            layoutReads.sortWith(
                Collections.reverseOrder(
                Comparator.comparingInt({ r: ReadLayout.Read -> r.alignedPosition + r.sequence.length })
                    .thenComparingDouble({ r: ReadLayout.Read -> r.baseQualities.average() }) // handle the highest quality ones first
                    .thenComparingInt({ r: ReadLayout.Read -> r.alignedPosition })
                    .thenComparing({ r: ReadLayout.Read -> r.readKey.readName }) // lastly we use read Id just in case
            ))
        }
        else
        {
            // start with sequences that are longest, highest quality first
            layoutReads.sortWith(
                Collections.reverseOrder(
                Comparator.comparingInt({ r: ReadLayout.Read -> r.alignedPosition })
                    .thenComparingDouble({ r: ReadLayout.Read -> r.baseQualities.average() }) // handle the highest quality ones first
                    .thenComparingInt({ r: ReadLayout.Read -> r.alignedPosition + r.sequence.length })
                    .thenComparing({ r: ReadLayout.Read -> r.readKey.readName }) // lastly we use read Id just in case
            ))
        }

        sLogger.info("building layouts from {} reads", layoutReads.size)

        val maxToCompareBeforeAlignedPos: Int = if (geneType.vj == VJ.V) CiderConstants.ASSUMED_ANCHOR_BASE_LENGTH else -1
        val maxToCompareAfterAlignedPos: Int = if (geneType.vj == VJ.J) CiderConstants.ASSUMED_ANCHOR_BASE_LENGTH else -1

        val readLayouts: MutableList<ReadLayout> = ArrayList()

        // go through the read data list, and add one by one to the list of clusters
        // if there are multiple clusters that matches, we choose the highest one
        for (readData in layoutReads)
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
                if (readMatchesLayout(readData, readLayout, minMatchedBases, minBaseQuality.toByte(),
                        maxToCompareBeforeAlignedPos, maxToCompareAfterAlignedPos))
                {
                    //sLogger.info("match found")
                    matchedLayout = readLayout
                }
            }

            if (matchedLayout == null)
            {
                val layout = ReadLayout()
                layout.addRead(readData, minBaseQuality.toByte())
                readLayouts.add(layout)
            }
            else
            {
                // add to all layouts that match
                // for (l in matchedLayouts)
                //addToOverlay(l, readData, minBaseQuality)
                matchedLayout.addRead(readData, minBaseQuality.toByte())
            }
        }

        sLogger.info("built {} layouts from {} reads", readLayouts.size, layoutReads.size)

        return readLayouts
    }

    companion object
    {
        private val sLogger = LogManager.getLogger(VJReadLayoutBuilder::class.java)

        //
        fun readMatchesLayout(readData: ReadLayout.Read, layout: ReadLayout,
                              minMatchedBases: Int, minBaseQuality: Byte,
                              maxToCompareBeforeAlignedPos: Int, // -1 means compare all
                              maxToCompareAfterAlignedPos: Int // -1 means compare all
        ): Boolean
        {
            require(maxToCompareBeforeAlignedPos >= 0 || maxToCompareBeforeAlignedPos == -1)
            require(maxToCompareAfterAlignedPos >= 0 || maxToCompareAfterAlignedPos == -1)

            // NOTE: alignedPosition is with respect to itself.
            // i.e. for a layout it is the position from layout start
            // for a read it is position from the read sequence start
            // it can actually be negative or larger than the sequence

            var startPosAdj = Math.min(layout.alignedPosition, readData.alignedPosition)

            if (maxToCompareBeforeAlignedPos != -1)
            {
                startPosAdj = Math.min(startPosAdj, maxToCompareBeforeAlignedPos)
            }

            // this will work also if the position is negative, i.e. before the sequence
            val layoutOffsetStart: Int = layout.alignedPosition - startPosAdj
            val readOffsetStart: Int = readData.alignedPosition - startPosAdj

            var layoutEnd: Int = layout.length
            var readEnd: Int = readData.readLength

            if (maxToCompareAfterAlignedPos != -1)
            {
                layoutEnd = Math.min(layoutEnd, layout.alignedPosition + maxToCompareAfterAlignedPos)
                readEnd = Math.min(readEnd, readData.alignedPosition + maxToCompareAfterAlignedPos)
            }

            // n is the number of bases overlap
            val n: Int = Math.min(layoutEnd - layoutOffsetStart, readEnd - readOffsetStart)

            if (n < minMatchedBases)
                return false

            var i1: Int = layoutOffsetStart
            var i2: Int = readOffsetStart
            var comparedCount: Int = 0
            var matchCount: Int = 0

            if (layoutOffsetStart + n > layout.length)
            {
                throw RuntimeException("read1StartOffset + n > read1Seq.length")
            }

            if (readOffsetStart + n > readData.readLength)
            {
                throw RuntimeException("read2StartOffset + n > read2Seq.length")
            }

            for (i in 0 until n)
            {
                val b1: Char = layout.highQualSequence[i1]
                val b2: Char = readData.sequence[i2]

                if (b1 == b2)
                {
                    // TODO: check with Peter
                    ++comparedCount
                    ++matchCount
                }
                else if ((readData.baseQualities[i2] >= minBaseQuality) &&
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
    }
}