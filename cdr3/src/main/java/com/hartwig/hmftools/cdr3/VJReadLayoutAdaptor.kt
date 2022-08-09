package com.hartwig.hmftools.cdr3

import com.hartwig.hmftools.cdr3.layout.ReadLayout
import com.hartwig.hmftools.cdr3.layout.ReadLayoutBuilder
import org.apache.logging.log4j.LogManager

// helper class to convert from the outer VJ classes to the layout classes
// create an interface to make it easier to test
interface IVJReadLayoutAdaptor
{
    fun getAnchorMatchType(layout: ReadLayout) : VJReadCandidate.AnchorMatchMethod
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
class VJReadLayoutAdaptor : IVJReadLayoutAdaptor
{
    private val sLogger = LogManager.getLogger(javaClass)

    fun readCandidateToLayoutRead(readCandidate: VJReadCandidate, trimBases: Int) : ReadLayout.Read
    {
        val retainStart: Int
        val retainEnd: Int
        val alignedPosition: Int

        if (readCandidate.vjGeneType.vj == VJ.V)
        {
            // for V, we layout from the anchor start, left to right
            // we are only interested in what comes after anchor start
            retainStart = Math.max(readCandidate.anchorOffsetStart, trimBases)
            retainEnd = readCandidate.readLength - trimBases

            // aligned position we must take into account that we remove all bases before vAnchor start
            alignedPosition = readCandidate.anchorOffsetEnd - 1 - retainStart
        }
        else
        {
            // for J, we layout from the anchor last, right to left
            // we are only interested in what comes before anchor end
            retainStart = trimBases
            retainEnd = Math.min(readCandidate.anchorOffsetEnd, readCandidate.readLength - trimBases)

            // aligned position we must take into account that we remove all bases before vAnchor start
            alignedPosition = readCandidate.anchorOffsetStart - retainStart
        }

        return ReadLayout.Read(readCandidate, ReadKey(readCandidate.read.readName, readCandidate.read.firstOfPairFlag),
            readCandidate.readSequence.substring(retainStart, retainEnd),
            readCandidate.baseQualities.sliceArray(retainStart until retainEnd),
            alignedPosition)
    }

    fun toReadCandidate(read: ReadLayout.Read) : VJReadCandidate
    {
        return read.source as VJReadCandidate
    }

    fun getReadCandidates(layout: ReadLayout) : List<VJReadCandidate>
    {
        return layout.reads.map({ read: ReadLayout.Read -> toReadCandidate(read)})
    }

    override fun getAnchorMatchType(layout: ReadLayout): VJReadCandidate.AnchorMatchMethod
    {
        return getReadCandidates(layout).first().anchorMatchMethod
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
            null

        if (anchorRange == null)
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
        return layout.consensusSequence().substring(getAnchorRange(geneType, layout)!!)
    }

    fun getAnchorSupport(geneType: VJGeneType, layout: ReadLayout) : String
    {
        return layout.highQualSupportString().substring(getAnchorRange(geneType, layout)!!)
    }

    fun getCdr3Sequence(geneType: VJGeneType, layout: ReadLayout) : String
    {
        return layout.consensusSequence().substring(getCdr3Range(geneType, layout)!!)
    }

    fun getCdr3Support(geneType: VJGeneType, layout: ReadLayout) : String
    {
        return layout.highQualSupportString().substring(getCdr3Range(geneType, layout)!!)
    }

    fun buildLayouts(geneType: VJGeneType, readCandidates: List<VJReadCandidate>,
                     minBaseQuality: Int, minMatchedBases: Int, minMatchRatio: Double,
                     numBasesToTrim: Int = 0)
    : List<ReadLayout>
    {
        sLogger.info("building {} layouts from {} reads", geneType, readCandidates.size)

        val layoutInputs = readCandidates.map({ o -> readCandidateToLayoutRead(o, numBasesToTrim) }).toList()
        // for V type we build from left, J we build from right
        val layoutBuilder = ReadLayoutBuilder(
            inputReads = layoutInputs,
            minBaseQuality = minBaseQuality,
            minMatchedBases = minMatchedBases,
            minMatchRatio = minMatchRatio,
            alignLeft = geneType.vj == VJ.V)
        val readOverlays = layoutBuilder.build()
        return readOverlays
    }
}