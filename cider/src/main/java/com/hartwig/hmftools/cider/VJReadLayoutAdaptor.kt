package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.layout.ReadLayout
import com.hartwig.hmftools.cider.layout.ReadLayoutBuilder
import htsjdk.samtools.SAMRecord
import org.apache.logging.log4j.LogManager

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
class VJReadLayoutAdaptor(private val trimBases: Int) : IVJReadLayoutAdaptor
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

    companion object
    {
        private val sLogger = LogManager.getLogger(VJReadLayoutAdaptor::class.java)

        private fun numTrailingPolyG(seq: String, sliceEnd: Int) : Int
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

        private fun numLeadingPolyC(seq: String, sliceStart: Int) : Int
        {
            for (i in sliceStart until seq.length)
            {
                if (seq[i] != 'C')
                {
                    return i
                }
            }
            return seq.length - sliceStart
        }
    }

    fun readCandidateToLayoutRead(readCandidate: VJReadCandidate) : ReadLayout.Read?
    {
        val slice = determineReadSlice(readCandidate.read, readCandidate.useReverseComplement,
                                        readCandidate.vjGeneType.vj, readCandidate.anchorOffsetStart,
                                        readCandidate.anchorOffsetEnd)

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
    private fun determineReadSlice(read: SAMRecord, useReverseComplement: Boolean, vj: VJ,
                            anchorOffsetStart: Int, anchorOffsetEnd: Int) : ReadSlice?
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

    // we also make sure they are aligned to codons
    fun getCdr3Range(geneType: VJGeneType, layout: ReadLayout) : IntRange?
    {
        val anchorRange = getAnchorRange(geneType, layout)

        if (anchorRange == null)
            return null

        // for V it is the sequence after the aligned position
        // for J it is the sequence before
        if (geneType.vj == VJ.V)
            return anchorRange.last + 1 until layout.consensusSequence().length
        else if (geneType.vj == VJ.J)
            return 0 until anchorRange.first
        return null
    }

    fun getAnchorSequence(geneType: VJGeneType, layout: ReadLayout) : String
    {
        val range = getAnchorRange(geneType, layout)
        return if (range != null)
            layout.consensusSequence().substring(range)
        else String()
    }

    fun getAnchorSupport(geneType: VJGeneType, layout: ReadLayout) : String
    {
        val range = getAnchorRange(geneType, layout)
        return if (range != null)
            layout.highQualSupportString().substring(range)
        else String()
    }

    fun getCdr3Sequence(geneType: VJGeneType, layout: ReadLayout) : String
    {
        val range = getCdr3Range(geneType, layout)
        return if (range != null)
            layout.consensusSequence().substring(range)
        else String()
    }

    fun getCdr3Support(geneType: VJGeneType, layout: ReadLayout) : String
    {
        val range = getCdr3Range(geneType, layout)
        return if (range != null)
            layout.highQualSupportString().substring(range)
        else String()
    }

    fun buildLayouts(geneType: VJGeneType, readCandidates: List<VJReadCandidate>,
                     minBaseQuality: Int, minMatchedBases: Int, minMatchRatio: Double)
    : List<ReadLayout>
    {
        sLogger.info("building {} layouts from {} reads", geneType, readCandidates.size)

        val layoutInputs = ArrayList<ReadLayout.Read>()

        for (r in readCandidates)
        {
            val layoutRead = readCandidateToLayoutRead(r)
            if (layoutRead != null)
                layoutInputs.add(layoutRead)
        }

        val layoutBuilder = ReadLayoutBuilder(
            inputReads = layoutInputs,
            minBaseQuality = minBaseQuality,
            minMatchedBases = minMatchedBases,
            minMatchRatio = minMatchRatio,
            alignLeft = geneType.vj == VJ.V         // for V type we build from left, J we build from right
        )
        val readLayouts = layoutBuilder.build()
        return readLayouts
    }
}