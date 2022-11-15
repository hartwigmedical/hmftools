package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.layout.ReadLayout
import org.apache.logging.log4j.LogManager

object VdjBuilderUtils
{
    private val sLogger = LogManager.getLogger(javaClass)

    data class SequenceOverlap(
        val seq1Offset: Int, // offset in overlap
        val seq2Offset: Int, // offset in overlap
        val overlapBases: Int,
        val highQualMatchBases: Int)

    // we try to overlap them. Note that it is possible for the V and J layout to overlap in any
    // direction
    //
    // Normally we have V layout on the left side
    // ----------------++++++++++++++++++++++         v layout
    //                 ++++++++++++++++++++++------   j layout
    //                 |--- overlap (>0) ---|
    //
    // it is also possible that
    //                 ++++++++++++++++++++++------   v layout
    // ----------------++++++++++++++++++++++         j layout
    //                 |--- overlap (<0) ---|
    //
    fun findSequenceOverlap(seq1: String, seq2: String, minOverlappedBases: Int) : SequenceOverlap?
    {
        var highQualMatchBases: Int = 0

        //  seq1            ==============
        //  seq2   ==================
        //  i      |--------|
        //  i is seq1Start - seq2Start
        //  i is negative in following:
        //  seq1   ==================
        //  seq2            ==================
        //  -i     |--------|
        for (i in -(seq1.length - minOverlappedBases) until (seq2.length - minOverlappedBases))
        {
            var seqMatch = true
            highQualMatchBases = 0
            val overlapSize = if (i > 0) Math.min(seq2.length - i, seq1.length) else Math.min(seq1.length + i, seq2.length)

            // check for overlap
            for (j in 0 until overlapSize)
            {
                val seq1Index: Int
                val seq2Index: Int

                if (i > 0)
                {
                    seq1Index = j
                    seq2Index = i + j
                }
                else
                {
                    seq1Index = -i + j
                    seq2Index = j
                }

                val vBase = seq1[seq1Index]
                val jBase = seq2[seq2Index]

                val bothHighQual = (vBase != 'N' && jBase != 'N')

                if (bothHighQual)
                {
                    if (vBase == jBase)
                    {
                        ++highQualMatchBases
                    }
                    else
                    {
                        seqMatch = false
                        break
                    }
                }
            }

            if (seqMatch && highQualMatchBases > minOverlappedBases)
            {
                // found overlap
                return if (i > 0)
                    SequenceOverlap(i, 0, overlapSize, highQualMatchBases)
                else
                    SequenceOverlap(0, -i, overlapSize, highQualMatchBases)
            }
        }

        return null
    }

    // we already worked out where the overlaps are, we now want to work out the
    // distance between V and J
    fun calcVJDistance(vSeq: String, jSeq: String, vPosInVSeq: Int, jPosInJSeq: Int, overlap: Int) : Int
    {
        // v layout:   TGC-GAATACC-CACATCCTGA-G
        // j layout:            CC-CACATCCTGA-GAGTG-TCAGA
        //                 |_____|            |___|
        //            V anchor(align)    J anchor(align)

        // v layout:    GC-GAATACC-CACATCCTGA-GAGTG-TCAGA
        // j layout: GATGC-GAATACC-CACATCCTGA-GAGTG-T
        //                 |_____|            |___|
        //            V anchor(align)    J anchor(align)
        //                        |----------|

        return 0
    }

    // TODO write unit test
    fun vdjSequenceCompare(vdj1: VDJSequence, vdj2: VDJSequence,
                                   diffAccumulator: (Map.Entry<Char, Int>, Map.Entry<Char, Int>) -> Boolean) : Boolean
    {
        require((vdj1.vAnchor != null) == (vdj2.vAnchor != null))
        require((vdj1.jAnchor != null) == (vdj2.jAnchor != null))

        // the lengths between V and J anchor must be equal
        if (vdj1.vAnchor != null && vdj1.jAnchor != null &&
            vdj2.vAnchor != null && vdj2.jAnchor != null &&
            vdj1.jAnchor.anchorBoundary - vdj1.vAnchor.anchorBoundary !=
            vdj2.jAnchor.anchorBoundary - vdj2.vAnchor.anchorBoundary)
        {
            return false
        }

        val (range1, range2) = getCompareRanges(vdj1, vdj2)

        // now we work out the ranges we can compare them
        var i1 = range1.first
        var i2 = range2.first

        while (i1 <= range1.last && i2 <= range2.last)
        {
            val support1 = vdj1.getSupportAt(i1)
            val support2 = vdj2.getSupportAt(i2)

            if (support1.key != support2.key)
            {
                if (!diffAccumulator(support1, support2))
                    return false
            }
            ++i1
            ++i2
        }
        return true
    }

    // get the ranges that line up both sequences based on anchor position
    fun getCompareRanges(vdj1: VDJSequence, vdj2: VDJSequence) : Pair<IntRange, IntRange>
    {
        require((vdj1.vAnchor != null) == (vdj2.vAnchor != null))
        require((vdj1.jAnchor != null) == (vdj2.jAnchor != null))

        // at least one side has to be there
        require(vdj1.vAnchor != null || vdj1.jAnchor != null)

        val offset: Int

        if (vdj1.vAnchor != null && vdj2.vAnchor != null)
        {
            offset = vdj1.vAnchor.anchorBoundary - vdj2.vAnchor.anchorBoundary
        }
        else if (vdj1.jAnchor != null && vdj2.jAnchor != null)
        {
            offset = vdj1.jAnchor.anchorBoundary - vdj2.jAnchor.anchorBoundary
        }
        else
        {
            // shouldn't get here
            throw RuntimeException()
        }

        var vdj1Start: Int = 0
        var vdj2Start: Int = 0

        // now they are lined up, we work out there the compare start / end is
        if (offset >= 0)
        {
            // vdj1 is further to the right
            vdj1Start = offset
        }
        else
        {
            // vdj2 is further to the right
            vdj2Start = -offset
        }

        val n = Math.min(vdj1.length - vdj1Start, vdj2.length - vdj2Start)
        var vdj1End: Int = vdj1Start + n
        var vdj2End: Int = vdj2Start + n

        if (vdj1.vAnchor == null && vdj2.vAnchor == null)
        {
            // if there is no V anchor, then we only compare up to 45 bases before J anchor
            vdj1Start = Math.max(vdj1.jAnchor!!.anchorBoundary - CiderConstants.PARTIAL_VDJ_UNANCHORED_LENGTH_BASES, vdj1Start)
            vdj2Start = Math.max(vdj2.jAnchor!!.anchorBoundary - CiderConstants.PARTIAL_VDJ_UNANCHORED_LENGTH_BASES, vdj2Start)
        }
        else if (vdj1.jAnchor == null && vdj2.jAnchor == null)
        {
            // if there is no J anchor, we only compare up to 45 bases after V anchor
            vdj1End = Math.min(vdj1.vAnchor!!.anchorBoundary + CiderConstants.PARTIAL_VDJ_UNANCHORED_LENGTH_BASES, vdj1End)
            vdj2End = Math.min(vdj2.vAnchor!!.anchorBoundary + CiderConstants.PARTIAL_VDJ_UNANCHORED_LENGTH_BASES, vdj2End)
        }

        assert(vdj1Start >= 0)
        assert(vdj2Start >= 0)
        assert(vdj1End <= vdj1.length)
        assert(vdj2End <= vdj2.length)

        val range1 = vdj1Start until vdj1End
        val range2 = vdj2Start until vdj2End

        return Pair(range1, range2)
    }

    // TODO: this function seems quite messy
    fun mergeVDJs(vdj1: VDJSequence, vdj2: VDJSequence, minBaseQuality: Byte): VDJSequence
    {
        require((vdj1.vAnchor != null) == (vdj2.vAnchor != null))
        require((vdj1.jAnchor != null) == (vdj2.jAnchor != null))

        val vdjPrimary: VDJSequence
        val vdjSecondary: VDJSequence

        if (vdj1.vAnchor != null && vdj1.jAnchor != null)
        {
            if (vdj1.jAnchor.anchorBoundary - vdj1.vAnchor.anchorBoundary != vdj2.jAnchor!!.anchorBoundary - vdj2.vAnchor!!.anchorBoundary)
            {
                sLogger.error(
                    "cannot merge {}(v:{}, j:{}) and {}(v:{}, j:{}), different lengths between anchors",
                    vdj1.aminoAcidSequenceFormatted, vdj1.vAnchor.matchMethod, vdj1.jAnchor.matchMethod,
                    vdj2.aminoAcidSequenceFormatted, vdj2.vAnchor.matchMethod, vdj2.jAnchor.matchMethod)

                throw RuntimeException("cannot merge VDJs: different lengths between anchors")
            }
        }

        sLogger.debug("start merge: {}(v:{}, j:{}, aligned pos:{}, within layout: {}-{}, v: {}, j: {}) and " +
                "{}(v:{}, j:{}, aligned pos:{}, within layout: {}-{}, v: {}, j: {})",
            vdj1.aminoAcidSequenceFormatted, vdj1.vAnchor?.matchMethod, vdj1.jAnchor?.matchMethod, vdj1.layout.alignedPosition,
            vdj1.layoutSliceStart, vdj1.layoutSliceEnd, vdj1.vAnchor?.anchorBoundary, vdj1.jAnchor?.anchorBoundary,
            vdj2.aminoAcidSequenceFormatted, vdj2.vAnchor?.matchMethod, vdj2.jAnchor?.matchMethod, vdj2.layout.alignedPosition,
            vdj2.layoutSliceStart, vdj2.layoutSliceEnd, vdj2.vAnchor?.anchorBoundary, vdj2.jAnchor?.anchorBoundary)

        val vAnchorBoundary = Math.max(vdj1.vAnchor?.anchorBoundary ?: 0, vdj2.vAnchor?.anchorBoundary ?: 0)
        val jAnchorBoundary = Math.max(vdj1.jAnchor?.anchorBoundary ?: 0, vdj2.jAnchor?.anchorBoundary ?: 0)

        // merge the anchor objects together
        val vAnchor: VJAnchor? = mergeAnchors(vdj1.vAnchor, vdj2.vAnchor, vdj1.vAnchorLength, vdj2.vAnchorLength, vAnchorBoundary)
        val jAnchor: VJAnchor? = mergeAnchors(vdj1.jAnchor, vdj2.jAnchor, vdj1.jAnchorLength, vdj2.jAnchorLength, jAnchorBoundary)

        val jAnchorLength: Int = Math.max(vdj1.jAnchorLength, vdj2.jAnchorLength)

        if (vdj1.numReads > vdj2.numReads)
        {
            vdjPrimary = vdj1
            vdjSecondary = vdj2
        }
        else
        {
            vdjPrimary = vdj2
            vdjSecondary = vdj1
        }

        // we calculate the aligned position for both sequences with respect to an anchor boundary
        // if there is V anchor, we use it, otherwise use J anchor. We just need to find a reference point
        val primaryAlignedPos = vdjPrimary.layout.alignedPosition - vdjPrimary.layoutSliceStart -
                (vdjPrimary.vAnchor?.anchorBoundary ?: vdjPrimary.jAnchor!!.anchorBoundary)
        val secondaryAlignedPos = vdjSecondary.layout.alignedPosition - vdjSecondary.layoutSliceStart -
                (vdjSecondary.vAnchor?.anchorBoundary ?: vdjSecondary.jAnchor!!.anchorBoundary)

        // this would shift reads to overlap each other, and create merged layout
        val secondaryLayoutAlignShift = primaryAlignedPos - secondaryAlignedPos
        val mergedLayout = ReadLayout.merge(vdjPrimary.layout, vdjSecondary.layout, 0, secondaryLayoutAlignShift, minBaseQuality)
        mergedLayout.id = "${vdjPrimary.layout.id};${vdjSecondary.layout.id}"

        // now we have to work out where to place the layout, using the final aligned position
        val boundaryPosDiff: Int = if (vAnchor != null)
            vdjPrimary.vAnchor!!.anchorBoundary - vAnchor.anchorBoundary
        else
            vdjPrimary.jAnchor!!.anchorBoundary - jAnchor!!.anchorBoundary


        //      vdj start                         vdj end
        //        |                                  |
        //    ++++-------V-------------------J-------++++    primary layout
        // +++++++-------V-------------------J-------+++++++ merged layout
        //        |                                  |
        //
        val posStartWithinLayout: Int = vdjPrimary.layoutSliceStart +
                boundaryPosDiff + // adjust for boundary position difference
                mergedLayout.alignedPosition - vdjPrimary.layout.alignedPosition // adjust for align position difference
        val posEndWithinLayout: Int

        require(vAnchor != null || jAnchor != null)

        if (jAnchor == null)
        {
            val postVLength = Math.max(vdjPrimary.length - vdjPrimary.vAnchorBoundary!!, vdjSecondary.length - vdjSecondary.vAnchorBoundary!!)
            posEndWithinLayout = posStartWithinLayout + vAnchor!!.anchorBoundary + postVLength
        }
        else
        {
            posEndWithinLayout = posStartWithinLayout + jAnchor!!.anchorBoundary + jAnchorLength
        }

        sLogger.debug("primary start: {}, primary v anchor: {}, v anchor: {}, layout aligned pos: {}, primary aligned pos: {}",
            vdjPrimary.layoutSliceStart, vdjPrimary.vAnchor?.anchorBoundary, vAnchor?.anchorBoundary,
            mergedLayout.alignedPosition, vdjPrimary.layout.alignedPosition)

        sLogger.debug("layout: {}, aligned pos: {}, pos within layout: {}-{}, v: {}, j: {}",
            mergedLayout.consensusSequence(), mergedLayout.alignedPosition, posStartWithinLayout, posEndWithinLayout,
            vAnchor?.anchorBoundary, jAnchor?.anchorBoundary)

        /*assert(supports.size > jAnchor.anchorBoundary)

        if (supports.size < jAnchorBoundary)
            throw RuntimeException()*/

        val combinedVDJ = VDJSequence(mergedLayout, posStartWithinLayout, posEndWithinLayout, vAnchor, jAnchor)

        sLogger.debug("merge: {}(v:{}, j:{}, read count: {}) and {}(v:{}, j:{}, read count: {}) together, merged: {}(v:{}, j:{}), 2ndary aligned pos shift: {}",
            vdj1.aminoAcidSequenceFormatted, vdj1.vAnchor?.matchMethod, vdj1.jAnchor?.matchMethod, vdj1.numReads,
            vdj2.aminoAcidSequenceFormatted, vdj2.vAnchor?.matchMethod, vdj2.jAnchor?.matchMethod, vdj2.numReads,
            combinedVDJ.aminoAcidSequenceFormatted, combinedVDJ.vAnchor?.matchMethod, combinedVDJ.jAnchor?.matchMethod,
            secondaryLayoutAlignShift)

        return combinedVDJ
    }


    fun mergeAnchors(anchor1: VJAnchor?, anchor2: VJAnchor?,
                     anchor1Length: Int, anchor2Length: Int,
                     newAnchorBoundary: Int) : VJAnchor?
    {
        if (anchor1 == null)
        {
            require(anchor2 == null)
            return null
        }

        if (anchor1 is VJAnchorByReadMatch)
        {
            if (anchor2 is VJAnchorByReadMatch)
            {
                // both are by read candidate match
                // choose the anchor that has higher number of reads
                val anchorWithMoreReads: VJAnchorByReadMatch = maxOf(anchor1, anchor1, Comparator.comparingInt(VJAnchorByReadMatch::numReads))

                // sum the number of reads together
                return anchorWithMoreReads.copy(anchorBoundary = newAnchorBoundary, numReads = anchor1.numReads + anchor2.numReads)
            }
            // copy anchor 1
            return anchor1.copy(anchorBoundary = newAnchorBoundary)
        }
        else if (anchor2 is VJAnchorByReadMatch)
        {
            // anchor1 is not by Read match
            // copy anchor 2
            return anchor2.copy(anchorBoundary = newAnchorBoundary)
        }

        // neither are by read match, we choose the longer one
        val anchorTemplate: VJAnchor? = if (anchor1Length >= anchor2Length) anchor1 else anchor2

        if (anchorTemplate is VJAnchorByReadMatch)
        {
            // smart casted to VJAnchorByReadMatch
            return anchorTemplate.copy(anchorBoundary = newAnchorBoundary)
        }
        else if (anchorTemplate is VJAnchorByBlosum)
        {
            // smart casted to VJAnchorByBlosum
            return anchorTemplate.copy(anchorBoundary = newAnchorBoundary)
        }

        throw IllegalArgumentException("unknown VJAnchor subtype: ${anchorTemplate!!::class}")
    }
}