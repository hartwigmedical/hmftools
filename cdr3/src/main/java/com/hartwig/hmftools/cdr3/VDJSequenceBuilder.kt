package com.hartwig.hmftools.cdr3

import com.google.common.collect.ArrayListMultimap
import com.hartwig.hmftools.cdr3.layout.ReadLayout
import com.hartwig.hmftools.cdr3.layout.VDJCandidate
import com.hartwig.hmftools.common.codon.Codons
import org.apache.logging.log4j.LogManager

object VDJSequenceBuilder
{
    private val sLogger = LogManager.getLogger(javaClass)

    @JvmStatic
    fun buildVDJSequences(layoutMultimap: Map<VJGeneType, List<ReadLayout>>, vjGeneStore: VJGeneStore, minBaseQuality: Byte) : List<VDJSequence>
    {
        val anchorBlosumSearcher = AnchorBlosumSearcher(
            vjGeneStore,
            CiderConstants.MIN_PARTIAL_ANCHOR_AA_LENGTH,
            CiderConstants.MAX_BLOSUM_DIFF_PER_AA,
            CiderConstants.ANCHOR_SIMILARITY_SCORE_CONSTANT
        )

        val vdjList = ArrayList<VDJSequence>()

        var seqIndex = 0
        for ((geneType, layouts) in layoutMultimap)
        {
            for (l in layouts)
            {
                val vdj = build(geneType, l, anchorBlosumSearcher, seqIndex)

                if (vdj != null)
                    vdjList.add(vdj)

                ++seqIndex
            }
        }

        // sort them by support min
        vdjList.sortByDescending({ vdj -> vdj.supportMin })

        sLogger.info("found {} vdj sequences before merge", vdjList.size)

        val vTypeVdjMap: ArrayListMultimap<VJGeneType, VDJSequence> = ArrayListMultimap.create()

        for (vdj in vdjList)
        {
            vTypeVdjMap.put(vdj.vAnchor.geneType, vdj)
        }

        vdjList.clear()

        for (vType in vTypeVdjMap.keySet())
        {
            vdjList.addAll(mergeIdentical(vTypeVdjMap.get(vType), minBaseQuality))
        }

        // sort again
        vdjList.sortByDescending({ vdj -> vdj.supportMin })

        sLogger.info("{} vdj sequences after merge", vdjList.size)

        return vdjList
    }

    private fun build(layoutGeneType: VJGeneType, layout: ReadLayout,
                      anchorBlosumSearcher: AnchorBlosumSearcher, seqIndex: Int)
        : VDJSequence?
    {
        val targetAnchorType = when (layoutGeneType)
        {
            VJGeneType.IGHV -> VJGeneType.IGHJ
            VJGeneType.IGHJ -> VJGeneType.IGHV
            VJGeneType.IGKV -> VJGeneType.IGKJ
            VJGeneType.IGKJ -> VJGeneType.IGKV
            VJGeneType.IGLV -> VJGeneType.IGLJ
            VJGeneType.IGLJ -> VJGeneType.IGLV
            VJGeneType.TRAV -> VJGeneType.TRAJ
            VJGeneType.TRAJ -> VJGeneType.TRAV
            VJGeneType.TRBV -> VJGeneType.TRBJ
            VJGeneType.TRBJ -> VJGeneType.TRBV
            VJGeneType.TRDV -> VJGeneType.TRDJ
            VJGeneType.TRDJ -> VJGeneType.TRDV
            VJGeneType.TRGV -> VJGeneType.TRGJ
            VJGeneType.TRGJ -> VJGeneType.TRGV
        }

        // we want to use the indices to work where things are
        val layoutSeq: String = layout.consensusSequence()

        val layoutAnchorRange: IntRange = VJReadLayoutAdaptor.getAnchorRange(layoutGeneType, layout)!!

        val searchStart = if (targetAnchorType.vj == VJ.V) 0 else layoutAnchorRange.last + 1
        val searchEnd = if (targetAnchorType.vj == VJ.V) layoutAnchorRange.first else layoutSeq.length

        // find all the homolog sequences
        val anchorBlosumMatches: List<AnchorBlosumSearcher.AnchorBlosumMatch> =
            anchorBlosumSearcher.findAnchorHomolog(layoutSeq, targetAnchorType, searchStart, searchEnd)

        if (anchorBlosumMatches.isEmpty())
            return null

        val anchorBlosumMatch = anchorBlosumMatches.first()

        var vAnchorRange: IntRange
        var jAnchorRange: IntRange

        if (layoutGeneType.vj == VJ.V)
        {
            vAnchorRange = VJReadLayoutAdaptor.getAnchorRange(layoutGeneType, layout)!!
            jAnchorRange = anchorBlosumMatch.anchorStart until anchorBlosumMatch.anchorEnd
        }
        else
        {
            vAnchorRange = anchorBlosumMatch.anchorStart until anchorBlosumMatch.anchorEnd
            jAnchorRange = VJReadLayoutAdaptor.getAnchorRange(layoutGeneType, layout)!!
        }

        val layoutStart = vAnchorRange.first
        val layoutEnd = jAnchorRange.last + 1 // inclusive to exclusive

        // make sure we shift the anchor ranges as well now we trimmed the left
        vAnchorRange = vAnchorRange.first - layoutStart .. vAnchorRange.last - layoutStart
        jAnchorRange = jAnchorRange.first - layoutStart .. jAnchorRange.last - layoutStart

        // construct a VDJ sequence
        val vAnchor: VJAnchor
        val jAnchor: VJAnchor

        if (layoutGeneType.vj == VJ.V)
        {
            vAnchor = createVJAnchorByReadMatch(
                type = VJAnchor.Type.V,
                anchorBoundary = vAnchorRange.last + 1,
                readCandidates = VJReadLayoutAdaptor.getReadCandidates(layout)
            )
            jAnchor = VJAnchorByBlosum(
                type = VJAnchor.Type.J,
                anchorBoundary = jAnchorRange.first,
                templateAnchorSeq = anchorBlosumMatch.templateAnchorSeq,
                templateGenes = anchorBlosumMatch.templateGenes,
                similarityScore = anchorBlosumMatch.similarityScore
            )
        }
        else
        {
            vAnchor = VJAnchorByBlosum(
                type = VJAnchor.Type.V,
                anchorBoundary = vAnchorRange.last + 1,
                templateAnchorSeq = anchorBlosumMatch.templateAnchorSeq,
                templateGenes = anchorBlosumMatch.templateGenes,
                similarityScore = anchorBlosumMatch.similarityScore
            )

            jAnchor = createVJAnchorByReadMatch(
                type = VJAnchor.Type.J,
                anchorBoundary = jAnchorRange.first,
                readCandidates = VJReadLayoutAdaptor.getReadCandidates(layout)
            )
        }

        val vdj = VDJSequence(seqIndex.toString(), layout, layoutStart, layoutEnd, vAnchor, jAnchor)

        sLogger.debug("built VDJ sequence: {}-{}-{}",
            Codons.aminoAcidFromBases(vdj.vAnchorSequence),
            Codons.aminoAcidFromBases(vdj.dSequenceShort),
            Codons.aminoAcidFromBases(vdj.jAnchorSequence))

        return vdj
    }

    // we should probably say allow some more leeway if a position has very little support
    fun tryOverlapVJ(vLayout: ReadLayout, jLayout: ReadLayout, minOverlappedBases: Int): VDJCandidate?
    {
        val vSupport = vLayout.highQualSequenceSupport.support
        val jSupport = jLayout.highQualSequenceSupport.support

        var overlappedBases: Int = -1

        // try multiple location to align them
        for (i in minOverlappedBases until Math.min(vSupport.size, jSupport.size))
        {
            var seqMatch = true

            // check for overlap
            for (j in 0 until i)
            {
                val vIndex = vSupport.size - i + j
                val jIndex = j

                val vBase = vSupport[vIndex].likelyBase
                val jBase = jSupport[jIndex].likelyBase

                if (vBase != 'N' && jBase != 'N' && vBase != jBase)
                {
                    seqMatch = false
                    break
                }
            }

            if (seqMatch)
            {
                // match found! now put them together into a VDJ sequence
                overlappedBases = i
                break
            }
        }

        if (overlappedBases == -1)
            return null

        /*val fullSeq = vSeq.substring(0, alignment.seq1AlignStart) + overlapSeqBuilder.toString() + jSeq.substring(alignment.seq2AlignEnd)

        sLogger.debug("v and j overlap")
        sLogger.debug("v: {}", vSeq)
        sLogger.debug("j: {}", jSeq)
        sLogger.debug("overlap: v({}), j({})", vOverlapSeq, jOverlapSeq)
        sLogger.debug("full seq: {}", fullSeq)

        return VDJCandidate(vLayout, jLayout, alignment.seq1AlignStart, alignment.seq2AlignEnd, overlapSeqBuilder.toString())*/
        return null
    }

    // for VDJ sequences that are identical, we put them together
    fun mergeIdentical(vdjList: List<VDJSequence>, minBaseQuality: Byte) : List<VDJSequence>
    {
        // first parse we do it in a fast way to find the ones that are exactly the same
        val seqHashMap: ArrayListMultimap<String, VDJSequence> = ArrayListMultimap.create()

        vdjList.forEach({ vdj -> seqHashMap.put(vdj.sequence, vdj) })

        /*
        val outVDJList = ArrayList<VDJSequence>()

        for (seq in seqHashMap.keySet())
        {
            val identicalVdjList = seqHashMap.get(seq)

            if (identicalVdjList.size == 1)
            {
                outVDJList.add(identicalVdjList[0])
                continue
            }

            var combinedVDJ = collapseVDJs(identicalVdjList[0], identicalVdjList[1])

            for (i in 1 until identicalVdjList.size)
            {
                combinedVDJ = collapseVDJs(combinedVDJ, identicalVdjList[i])
            }
            outVDJList.add(combinedVDJ)
        }*/

        val outVDJList = ArrayList<VDJSequence>(vdjList)

        // now we again compare them all against one another. The reason is two
        // 1. the sequences could be subsequence of another
        // 2. we also want to remove sequences that are nearly identical
        var i: Int = 0
        var j: Int = 1
        while (i < outVDJList.size && j < outVDJList.size)
        {
            val vdj1 = outVDJList[i]
            val vdj2 = outVDJList[j]
            var numBaseDiff = 0
            val diffAccumulator = { s1: Map.Entry<Char, Int>, s2: Map.Entry<Char, Int> ->
                if (s1.value >= 1 && s2.value >= 1)
                {
                    ++numBaseDiff
                }
                numBaseDiff == 0
            }
            // count them as different if one base with one read is different
            if (vdjEqual(vdj1, vdj2, diffAccumulator))
            {
                val combinedVDJ = mergeVDJs(vdj1, vdj2, minBaseQuality)

                outVDJList[i] = combinedVDJ
                // remove vdj2 from the list
                outVDJList.removeAt(j)
                --j
            }

            ++j

            if (j >= outVDJList.size)
            {
                ++i
                j = i + 1
            }
        }

        // now compare all those that are just one base diff
        i = 0
        j = 1
        while (i < outVDJList.size && j < outVDJList.size)
        {
            val vdj1 = outVDJList[i]
            val vdj2 = outVDJList[j]

            var numBaseDiff = 0
            var maxReadsAtDiff = 0
            val diffAccumulator = { s1: Map.Entry<Char, Int>, s2: Map.Entry<Char, Int> ->
                if (s1.value >= 1 && s2.value >= 1)
                {
                    ++numBaseDiff
                    maxReadsAtDiff = Math.max(maxReadsAtDiff, Math.min(s1.value, s2.value))
                }

                // we say they are the same if there is only 1 base difference
                // and only 1 read supporting that one base diff
                (numBaseDiff <= 1 && maxReadsAtDiff < 2)
            }

            if (vdjEqual(vdj1, vdj2, diffAccumulator))
            {
                sLogger.debug("removing {} as it is similar to {}",
                    vdj2.aminoAcidSequenceFormatted, vdj1.aminoAcidSequenceFormatted)

                // remove vdj2 from the list
                outVDJList.removeAt(j)
                --j
            }

            ++j

            if (j >= outVDJList.size)
            {
                ++i
                j = i + 1
            }
        }

        return outVDJList
    }

    internal fun mergeVDJs(vdj1: VDJSequence, vdj2: VDJSequence, minBaseQuality: Byte): VDJSequence
    {
        val vdjPrimary: VDJSequence
        val vdjSecondary: VDJSequence

        if (vdj1.jAnchor.anchorBoundary - vdj1.vAnchor.anchorBoundary != vdj2.jAnchor.anchorBoundary - vdj2.vAnchor.anchorBoundary)
        {
            sLogger.error("cannot merge {}(v:{}, j:{}) and {}(v:{}, j:{}), different lengths between anchors",
                vdj1.aminoAcidSequenceFormatted, vdj1.vAnchor.matchMethod, vdj1.jAnchor.matchMethod,
                vdj2.aminoAcidSequenceFormatted, vdj2.vAnchor.matchMethod, vdj2.jAnchor.matchMethod)

            throw RuntimeException("cannot merge VDJs: different lengths between anchors")
        }

        sLogger.debug("start merge: {}(v:{}, j:{}, aligned pos:{}, within layout: {}-{}, v: {}, j: {}) and " +
                "{}(v:{}, j:{}, aligned pos:{}, within layout: {}-{}, v: {}, j: {})",
            vdj1.aminoAcidSequenceFormatted, vdj1.vAnchor.matchMethod, vdj1.jAnchor.matchMethod, vdj1.layout.alignedPosition,
            vdj1.layoutStart, vdj1.layoutEnd, vdj1.vAnchor.anchorBoundary, vdj1.jAnchor.anchorBoundary,
            vdj2.aminoAcidSequenceFormatted, vdj2.vAnchor.matchMethod, vdj2.jAnchor.matchMethod, vdj2.layout.alignedPosition,
            vdj2.layoutStart, vdj2.layoutEnd, vdj2.vAnchor.anchorBoundary, vdj2.jAnchor.anchorBoundary)

        if (vdj1.supportMin > vdj2.supportMin)
        {
            vdjPrimary = vdj1
            vdjSecondary = vdj2
        }
        else
        {
            vdjPrimary = vdj2
            vdjSecondary = vdj1
        }

        val vAnchorBoundary = Math.max(vdj1.vAnchor.anchorBoundary, vdj2.vAnchor.anchorBoundary)
        val jAnchorBoundary = Math.max(vdj1.jAnchor.anchorBoundary, vdj2.jAnchor.anchorBoundary)

        // we first prefer anchor that is longer, then we prefer anchor
        // that is by read match
        val vAnchorTemplate: VJAnchor
        val jAnchorTemplate: VJAnchor
        val jAnchorLength: Int

        if (vdjPrimary.vAnchorLength < vdjSecondary.vAnchorLength ||
            vdjPrimary.vAnchor !is VJAnchorByReadMatch && vdjSecondary.vAnchor is VJAnchorByReadMatch)
        {
            vAnchorTemplate = vdjSecondary.vAnchor
        }
        else
        {
            vAnchorTemplate = vdjPrimary.vAnchor
        }
        if (vdjPrimary.jAnchorLength < vdjSecondary.jAnchorLength ||
            vdjPrimary.jAnchor !is VJAnchorByReadMatch && vdjSecondary.jAnchor is VJAnchorByReadMatch)
        {
            jAnchorTemplate = vdjSecondary.jAnchor
            jAnchorLength = vdjSecondary.jAnchorLength
        }
        else
        {
            jAnchorTemplate = vdjPrimary.jAnchor
            jAnchorLength = vdjPrimary.jAnchorLength
        }

        // create the V and J anchors from template
        val vAnchor: VJAnchor = copyVjAnchorTemplate(vAnchorTemplate, vAnchorBoundary)
        val jAnchor: VJAnchor = copyVjAnchorTemplate(jAnchorTemplate, jAnchorBoundary)

        // we calculate the aligned position for both sequences with respect to the V anchor boundary
        val primaryAlignedPos = vdjPrimary.layout.alignedPosition - vdjPrimary.layoutStart - vdjPrimary.vAnchor.anchorBoundary
        val secondaryAlignedPos = vdjSecondary.layout.alignedPosition - vdjSecondary.layoutStart - vdjSecondary.vAnchor.anchorBoundary

        // this would shift them to overlap each other
        val secondaryLayoutOffsetShift = primaryAlignedPos - secondaryAlignedPos
        val layout = ReadLayout.merge(vdjPrimary.layout, vdjSecondary.layout, 0, secondaryLayoutOffsetShift, minBaseQuality)

        // now we have to work out where to place the layout, using the final aligned position
        val posStartWithinLayout = vdjPrimary.layoutStart +
                vdjPrimary.vAnchor.anchorBoundary - vAnchor.anchorBoundary + // adjust for V anchor length
                layout.alignedPosition - vdjPrimary.layout.alignedPosition // adjust for align position difference
        val posEndWithinLayout = posStartWithinLayout + jAnchor.anchorBoundary + jAnchorLength

        sLogger.debug("primary start: {}, primary v anchor: {}, v anchor: {}, layout aligned pos: {}, primary aligned pos: {}",
            vdjPrimary.layoutStart, vdjPrimary.vAnchor.anchorBoundary, vAnchor.anchorBoundary,
                    layout.alignedPosition, vdjPrimary.layout.alignedPosition)

        sLogger.debug("layout: {}, aligned pos: {}, pos within layout: {}-{}, v: {}, j: {}",
            layout.consensusSequence(), layout.alignedPosition, posStartWithinLayout, posEndWithinLayout,
            vAnchor.anchorBoundary, jAnchor.anchorBoundary)

        /*assert(supports.size > jAnchor.anchorBoundary)

        if (supports.size < jAnchorBoundary)
            throw RuntimeException()*/

        val combinedVDJ = VDJSequence(vdjPrimary.id, layout, posStartWithinLayout, posEndWithinLayout, vAnchor, jAnchor)

        sLogger.debug("merge: {}(v:{}, j:{}, read count: {}) and {}(v:{}, j:{}, read count: {}) together, merged: {}(v:{}, j:{}), 2ndary aligned pos shift: {}",
            vdj1.aminoAcidSequenceFormatted, vdj1.vAnchor.matchMethod, vdj1.jAnchor.matchMethod, vdj1.numReads,
            vdj2.aminoAcidSequenceFormatted, vdj2.vAnchor.matchMethod, vdj2.jAnchor.matchMethod, vdj2.numReads,
            combinedVDJ.aminoAcidSequenceFormatted, combinedVDJ.vAnchor.matchMethod, combinedVDJ.jAnchor.matchMethod,
            secondaryLayoutOffsetShift)

        return combinedVDJ
    }

    // TODO write unit test
    private fun vdjEqual(vdj1: VDJSequence, vdj2: VDJSequence,
                         diffAccumulator: (Map.Entry<Char, Int>, Map.Entry<Char, Int>) -> Boolean) : Boolean
    {
        // the lengths between V and J anchor must be equal
        if (vdj1.jAnchor.anchorBoundary - vdj1.vAnchor.anchorBoundary != vdj2.jAnchor.anchorBoundary - vdj2.vAnchor.anchorBoundary)
        {
            return false
        }

        val (range1, range2) = getVdjRanges(vdj1, vdj2)

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
    fun getVdjRanges(vdj1: VDJSequence, vdj2: VDJSequence) : Pair<IntRange, IntRange>
    {
        var vdj1Start: Int = 0
        var vdj1End: Int = vdj1.length
        var vdj2Start: Int = 0
        var vdj2End: Int = vdj2.length

        if (vdj1.vAnchor.anchorBoundary > vdj2.vAnchor.anchorBoundary)
        {
            vdj1Start = vdj1.vAnchor.anchorBoundary - vdj2.vAnchor.anchorBoundary
        }
        else if (vdj1.vAnchor.anchorBoundary < vdj2.vAnchor.anchorBoundary)
        {
            vdj2Start = vdj2.vAnchor.anchorBoundary - vdj1.vAnchor.anchorBoundary
        }

        val endDiff = (vdj1.length - vdj1Start) - (vdj2.length - vdj2Start)

        if (endDiff > 0)
        {
            vdj1End -= endDiff
        }
        else if (endDiff < 0)
        {
            vdj2End -= endDiff
        }

        assert(vdj1Start >= 0)
        assert(vdj2Start >= 0)

        val range1 = vdj1Start until vdj1End
        val range2 = vdj2Start until vdj2End

        return Pair(range1, range2)
    }

    fun copyVjAnchorTemplate(anchorTemlpate: VJAnchor, anchorBoundary: Int): VJAnchor
    {
        if (anchorTemlpate is VJAnchorByReadMatch)
        {
            return anchorTemlpate.copy(anchorBoundary = anchorBoundary)
        }
        else if (anchorTemlpate is VJAnchorByBlosum)
        {
            return anchorTemlpate.copy(anchorBoundary = anchorBoundary)
        }
        else
        {
            throw RuntimeException("unknown VJAnchor child type")
        }
    }

    fun createVJAnchorByReadMatch(
        type: VJAnchor.Type,
        anchorBoundary: Int,
        readCandidates: List<VJReadCandidate>) : VJAnchorByReadMatch
    {
        val matchMethod =
            when (readCandidates.first().anchorMatchType)
            {
                VJReadCandidate.AnchorMatchType.ALIGN -> "align"
                VJReadCandidate.AnchorMatchType.EXACT -> "exact"
            }

        val geneTypes: List<VJGeneType> = readCandidates.map({ o -> o.vjGeneType }).distinct()
        if (geneTypes.size != 1)
            throw IllegalStateException("VJAnchorByReadMatch: gene types(${geneTypes}) size != 1")
        val geneType = geneTypes.first()

        return VJAnchorByReadMatch(type, geneType, anchorBoundary, matchMethod)
    }
}