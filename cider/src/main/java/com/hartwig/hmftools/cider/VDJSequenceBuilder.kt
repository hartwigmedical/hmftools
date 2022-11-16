package com.hartwig.hmftools.cider

import com.google.common.collect.ArrayListMultimap
import com.hartwig.hmftools.cider.layout.ReadLayout
import com.hartwig.hmftools.common.codon.Codons
import org.apache.logging.log4j.LogManager

class VDJSequenceBuilder(private val vjLayoutAdaptor: IVJReadLayoutAdaptor,
                         private val anchorBlosumSearcher: IAnchorBlosumSearcher,
                         private val minBaseQuality: Byte,
                         private val minOverlappedBases: Int)
{
    fun buildVDJSequences(layoutMultimap: Map<VJGeneType, List<ReadLayout>>) : List<VDJSequence>
    {
        val vdjList = ArrayList<VDJSequence>()

        // remember the layouts which are not used
        val oneSidedLayouts: ArrayListMultimap<VJGeneType, ReadLayout> = ArrayListMultimap.create()

        for ((geneType, layouts) in layoutMultimap)
        {
            for (l in layouts)
            {
                oneSidedLayouts.put(geneType, l)
            }
        }

        // try to join together the layouts
        for (vGeneType: VJGeneType in oneSidedLayouts.keySet().filter({ vjGeneType -> vjGeneType.vj == VJ.V }))
        {
            val vLayoutItr: MutableIterator<ReadLayout> = oneSidedLayouts[vGeneType].iterator()

            // we try each pair see if they can be aligned together
            // we find the best aligned pair and remove them and try again
            while (vLayoutItr.hasNext())
            {
                val vLayout: ReadLayout = vLayoutItr.next()
                loop@ for (jGeneType: VJGeneType in vGeneType.pairedVjGeneTypes())
                {
                    // try all the layouts that are paired types
                    val jLayoutItr: MutableIterator<ReadLayout> = oneSidedLayouts[jGeneType].iterator()
                    while (jLayoutItr.hasNext())
                    {
                        val jLayout: ReadLayout = jLayoutItr.next()
                        val vdj: VDJSequence? = tryOverlapVJ(vLayout, jLayout, vGeneType, jGeneType)

                        if (vdj != null)
                        {
                            vdjList.add(vdj)
                            vLayoutItr.remove()
                            jLayoutItr.remove()
                            break@loop // break out of the outer loop
                        }
                    }
                }
            }
        }

        // for the remaining ones, we try to complete with blosum
        val itr = oneSidedLayouts.entries().iterator()
        while (itr.hasNext())
        {
            val (geneType, layout) = itr.next()
            val vdj = tryCompleteLayoutWithBlosum(geneType, layout)

            if (vdj != null)
            {
                vdjList.add(vdj)
                itr.remove()
            }
        }

        // add all of the VDJ with only one side
        for ((geneType, layout) in oneSidedLayouts.entries())
        {
            val vdj = tryCreateOneSidedVdj(geneType, layout)
            if (vdj != null)
            {
                vdjList.add(vdj)
            }
        }

        // sort them by number of reads
        vdjList.sortByDescending({ vdj -> vdj.numReads })

        sLogger.info("found {} vdj sequences before merge", vdjList.size)

        val vjGeneToVdjMap: ArrayListMultimap<Pair<VJGeneType?, VJGeneType?>, VDJSequence> = ArrayListMultimap.create()

        for (vdj in vdjList)
        {
            if (vdj.vAnchor == null && vdj.jAnchor == null)
            {
                assert(false)
                continue
            }
            vjGeneToVdjMap.put(Pair(vdj.vAnchor?.geneType, vdj.jAnchor?.geneType), vdj)
        }

        vdjList.clear()

        for (vjType in vjGeneToVdjMap.keySet())
        {
            vdjList.addAll(mergeIdentical(vjGeneToVdjMap.get(vjType), minBaseQuality))
        }

        // sort again
        vdjList.sortByDescending({ vdj -> vdj.numReads })

        sLogger.info("{} vdj sequences after merge", vdjList.size)

        return vdjList
    }

    // This function tries to complete the layout by searching for the missing anchor
    // we already have one side of the anchor
    // use the blosum searcher to find the other side
    // if it is found then build a VDJSequence
    fun tryCompleteLayoutWithBlosum(layoutGeneType: VJGeneType, layout: ReadLayout)
            : VDJSequence?
    {
        sLogger.debug("try complete {} layout: {}", layoutGeneType, layout.consensusSequence())

        val targetAnchorTypes = layoutGeneType.pairedVjGeneTypes()

        // we want to use the indices to work where things are
        val layoutSeq: String = layout.consensusSequence()

        val layoutAnchorRange: IntRange? = vjLayoutAdaptor.getAnchorRange(layoutGeneType.vj, layout)

        if (layoutAnchorRange == null)
            return null

        val searchStart: Int
        val searchEnd: Int

        if (layoutGeneType.vj == VJ.V)
        {
            searchStart = layoutAnchorRange.last + 1
            searchEnd = layoutSeq.length
        }
        else
        {
            searchStart = 0
            searchEnd = layoutAnchorRange.first
        }

        // find all the homolog sequences
        val anchorBlosumMatch: AnchorBlosumMatch? = anchorBlosumSearcher.searchForAnchor(
                layoutSeq, targetAnchorTypes,
                IAnchorBlosumSearcher.Mode.DISALLOW_NEG_SIMILARITY,
                searchStart, searchEnd)

        if (anchorBlosumMatch == null)
            return null

        var vAnchorRange: IntRange?
        var jAnchorRange: IntRange?

        if (layoutGeneType.vj == VJ.V)
        {
            vAnchorRange = vjLayoutAdaptor.getAnchorRange(layoutGeneType.vj, layout)
            jAnchorRange = anchorBlosumMatch.anchorStart until anchorBlosumMatch.anchorEnd
        }
        else
        {
            vAnchorRange = anchorBlosumMatch.anchorStart until anchorBlosumMatch.anchorEnd
            jAnchorRange = vjLayoutAdaptor.getAnchorRange(layoutGeneType.vj, layout)
        }

        if (vAnchorRange == null || jAnchorRange == null)
            return null

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
                anchorBoundary = vAnchorRange.last + 1,
                vjGeneType = layoutGeneType,
                layout = layout
            )
            jAnchor = VJAnchorByBlosum(
                vj = VJ.J,
                geneType = anchorBlosumMatch.templateGenes.first().type,
                anchorBoundary = jAnchorRange.first,
                templateAnchorSeq = anchorBlosumMatch.templateAnchorSeq,
                templateGenes = anchorBlosumMatch.templateGenes,
                similarityScore = anchorBlosumMatch.similarityScore
            )
        }
        else
        {
            vAnchor = VJAnchorByBlosum(
                vj = VJ.V,
                geneType = anchorBlosumMatch.templateGenes.first().type,
                anchorBoundary = vAnchorRange.last + 1,
                templateAnchorSeq = anchorBlosumMatch.templateAnchorSeq,
                templateGenes = anchorBlosumMatch.templateGenes,
                similarityScore = anchorBlosumMatch.similarityScore
            )

            jAnchor = createVJAnchorByReadMatch(
                anchorBoundary = jAnchorRange.first,
                vjGeneType = layoutGeneType,
                layout = layout
            )
        }

        val vdj = VDJSequence(layout, layoutStart, layoutEnd, vAnchor, jAnchor)

        sLogger.debug("built VDJ sequence: {}, {}",
            vdj.aminoAcidSequenceFormatted, vdj.sequenceFormatted)

        return vdj
    }

    // we should probably say allow some more leeway if a position has very little support
    fun tryOverlapVJ(vLayout: ReadLayout, jLayout: ReadLayout, vLayoutGeneType: VJGeneType, jLayoutGeneType: VJGeneType): VDJSequence?
    {
        assert(vLayoutGeneType.vj == VJ.V)
        assert(jLayoutGeneType.vj == VJ.J)

        if (!vLayoutGeneType.pairedVjGeneTypes().contains(jLayoutGeneType))
        {
            // mismatched type, i.e. IGHV must pair with IGHJ
            return null
        }

        val vLayoutSeq = vLayout.highQualSequence
        val jLayoutSeq = jLayout.highQualSequence

        val vjLayoutOverlap: VdjBuilderUtils.SequenceOverlap? = VdjBuilderUtils.findSequenceOverlap(vLayoutSeq, jLayoutSeq, minOverlappedBases)

        if (vjLayoutOverlap == null)
            return null

        // consider the following
        //                  v aligned pos    j aligned pos
        //                      |                |
        //  ---------------VVVVV-------------               v layout
        //                  +++++++++++++++++----JJJJJ----  j layout
        //  ---------------VVVVV-----------------JJJJJ----  combine layout
        //  <--------------><--------------->
        //   jLayoutOffset       overlap
        //                      <----------------
        //                     jAlignedPositionShift
        //
        // work out how much j aligned position needs to be shifted so that they are all aligned in the
        // v aligned position
        val jAlignedPositionShift = (vjLayoutOverlap.seq1Offset + vLayout.alignedPosition) -
                                    (vjLayoutOverlap.seq2Offset + jLayout.alignedPosition)

        val combinedVjLayout: ReadLayout = ReadLayout.merge(jLayout, vLayout, jAlignedPositionShift, 0, minBaseQuality)
        combinedVjLayout.id = "${vLayout.id},${jLayout.id}"

        // next we want to work out where the anchor ranges are in the combined layout
        var vAnchorRange: IntRange? = vjLayoutAdaptor.getAnchorRange(VJ.V, vLayout)
        var jAnchorRange: IntRange? = vjLayoutAdaptor.getAnchorRange(VJ.J, jLayout)

        if (vAnchorRange == null || jAnchorRange == null)
            return null

        // now we use the aligned position of the combined layout
        val vAnchorShift = combinedVjLayout.alignedPosition - vLayout.alignedPosition
        vAnchorRange = vAnchorRange.first + vAnchorShift .. vAnchorRange.last + vAnchorShift

        val jAnchorShift = combinedVjLayout.alignedPosition - jLayout.alignedPosition - jAlignedPositionShift
        jAnchorRange = jAnchorRange.first + jAnchorShift .. jAnchorRange.last + jAnchorShift

        sLogger.debug("overlap: {}, jAlignShift: {}, vLayout: {}, jLayout: {}, combinedLayout: {}",
            vjLayoutOverlap, jAlignedPositionShift,
            vLayout.consensusSequence(), jLayout.consensusSequence(), combinedVjLayout.consensusSequence())

        val layoutSliceStart: Int
        val layoutSliceEnd = jAnchorRange.last + 1 // inclusive to exclusive

        var vAnchorBoundary = vAnchorRange.last + 1
        var jAnchorBoundary = jAnchorRange.first

        // construct a VDJ sequence
        val vAnchor: VJAnchor?
        val jAnchor: VJAnchor

        //if (jAnchorBoundary - vAnchorBoundary >= 0)
        if (vAnchorRange.first < jAnchorRange.first &&
            vAnchorRange.last < jAnchorRange.last)
        {
            layoutSliceStart = vAnchorRange.first
            vAnchorBoundary -= layoutSliceStart
            jAnchorBoundary -= layoutSliceStart

            vAnchor = createVJAnchorByReadMatch(
                anchorBoundary = vAnchorBoundary,
                vjGeneType = vLayoutGeneType,
                layout = vLayout)
        }
        else
        {
            // if the V and J anchor overlap, we remove the V anchor
            // reason is that J rearrangement happens more. But maybe we should rethink this

            vAnchor = null
            // we no longer trim to the v anchor start
            layoutSliceStart = 0
            sLogger.debug("VDJ vAnchorBoundary({}) < jAnchorBoundary({}), setting v anchor to null",
                vAnchorBoundary, jAnchorBoundary)
        }

        jAnchor = createVJAnchorByReadMatch(
            anchorBoundary = jAnchorBoundary,
            vjGeneType = jLayoutGeneType,
            layout = jLayout)

        val vdj = VDJSequence(combinedVjLayout, layoutSliceStart, layoutSliceEnd, vAnchor, jAnchor)

        sLogger.debug("built VDJ sequence: {}, by overlapping V layout({}): {} and J layout({}): {}",
            vdj.aminoAcidSequenceFormatted,
            vLayout.id, vLayout.consensusSequence(), jLayout.id, jLayout.consensusSequence())

        return vdj
    }

    fun tryCreateOneSidedVdj(layoutGeneType: VJGeneType, layout: ReadLayout): VDJSequence?
    {
        sLogger.debug("create one sided {} layout: {}", layoutGeneType, layout.consensusSequence())

        // we want to use the indices to work where things are
        val layoutAnchorRange: IntRange? = vjLayoutAdaptor.getAnchorRange(layoutGeneType.vj, layout)

        if (layoutAnchorRange == null)
            return null

        // construct a VDJ sequence
        val vAnchor: VJAnchor?
        val jAnchor: VJAnchor?
        val vdj: VDJSequence

        if (layoutGeneType.vj == VJ.V)
        {
            val vAnchorRange = vjLayoutAdaptor.getAnchorRange(layoutGeneType.vj, layout)

            if (vAnchorRange == null)
                return null

            val layoutStart = vAnchorRange.first

            // we limit one sided VDJ post V anchor length
            val layoutEnd = Math.min(layout.length, vAnchorRange.last + 1 + CiderConstants.PARTIAL_VDJ_UNANCHORED_LENGTH_BASES)

            vAnchor = createVJAnchorByReadMatch(
                anchorBoundary = vAnchorRange.last + 1 - layoutStart,
                vjGeneType = layoutGeneType,
                layout = layout
            )
            vdj = VDJSequence(layout, layoutStart, layoutEnd, vAnchor, null)

            sLogger.debug("built V only sequence: {}-{}",
                Codons.aminoAcidFromBases(vdj.vAnchorSequence),
                Codons.aminoAcidFromBases(vdj.cdr3SequenceShort))
        }
        else
        {
            val jAnchorRange = vjLayoutAdaptor.getAnchorRange(layoutGeneType.vj, layout)

            if (jAnchorRange == null)
                return null

            // we limit one sided VDJ pre J anchor length
            val layoutStart = Math.max(0, jAnchorRange.first - CiderConstants.PARTIAL_VDJ_UNANCHORED_LENGTH_BASES)
            val layoutEnd = jAnchorRange.last + 1

            jAnchor = createVJAnchorByReadMatch(
                anchorBoundary = jAnchorRange.first - layoutStart,
                vjGeneType = layoutGeneType,
                layout = layout
            )
            vdj = VDJSequence(layout, layoutStart, layoutEnd, null, jAnchor)

            sLogger.debug("built J only sequence: {}-{}",
                Codons.aminoAcidFromBases(vdj.cdr3SequenceShort),
                Codons.aminoAcidFromBases(vdj.jAnchorSequence))
        }

        return vdj
    }

    fun createVJAnchorByReadMatch(
        anchorBoundary: Int,
        layout: ReadLayout,
        vjGeneType: VJGeneType) : VJAnchorByReadMatch
    {
        val matchMethodStr =
            when (vjLayoutAdaptor.getAnchorMatchMethod(layout))
            {
                VJReadCandidate.MatchMethod.ALIGN -> "align"
                VJReadCandidate.MatchMethod.EXACT -> "exact"
                VJReadCandidate.MatchMethod.BLOSUM -> "blosum"
            }

        val templateAnchorSeq = vjLayoutAdaptor.getTemplateAnchorSequence(layout)

        return VJAnchorByReadMatch(vjGeneType.vj, vjGeneType, anchorBoundary, matchMethodStr,
            templateAnchorSeq, layout.reads.size)
    }

    companion object
    {
        private val sLogger = LogManager.getLogger(VDJSequenceBuilder::class.java)

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
                // count them as different if one base with one read is different
                if (vdjSequenceIdentical(vdj1, vdj2))
                {
                    val combinedVDJ = VdjBuilderUtils.mergeVDJs(vdj1, vdj2, minBaseQuality)

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

                if (VdjBuilderUtils.vdjSequenceCompare(vdj1, vdj2, diffAccumulator))
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

        fun vdjSequenceIdentical(vdj1: VDJSequence, vdj2: VDJSequence) : Boolean
        {
            var numBaseDiff = 0
            val diffAccumulator = { baseSupport1: Map.Entry<Char, Int>, baseSupport2: Map.Entry<Char, Int> ->

                // only count diff if support is not zero
                if (baseSupport1.value >= 1 && baseSupport2.value >= 1)
                {
                    ++numBaseDiff
                }

                // we keep comparing if there is no diff
                numBaseDiff == 0// return
            }
            return VdjBuilderUtils.vdjSequenceCompare(vdj1, vdj2, diffAccumulator)
        }
    }
}