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

        var seqIndex = 0
        for ((geneType, layouts) in layoutMultimap)
        {
            for (l in layouts)
            {
                val vdj = tryCompleteLayoutWithBlosum(geneType, l, seqIndex)

                if (vdj != null)
                {
                    vdjList.add(vdj)
                }
                else
                {
                    oneSidedLayouts.put(geneType, l)
                }

                ++seqIndex
            }
        }

        // try to join together the orphaned layouts
        for (vGeneType: VJGeneType in oneSidedLayouts.keySet().filter({ vjGeneType -> vjGeneType.vj == VJ.V }))
        {
            val jGeneType: VJGeneType = CiderUtils.getPairedVjGeneType(vGeneType)
            val vLayoutList: MutableList<ReadLayout> = oneSidedLayouts[vGeneType]
            val jLayoutList: MutableList<ReadLayout> = oneSidedLayouts[jGeneType]

            // we try each pair see if they can be aligned together
            // we find the best aligned pair and remove them and try again
            var i = 0
            while (i < vLayoutList.size)
            {
                val vLayout: ReadLayout = vLayoutList[i]
                var j = 0
                while (j < jLayoutList.size)
                {
                    val jLayout: ReadLayout = jLayoutList[j]
                    val vdj: VDJSequence? = tryOverlapVJ(vLayout, jLayout, vGeneType, jGeneType)

                    if (vdj != null)
                    {
                        vdjList.add(vdj)
                        vLayoutList.removeAt(i)
                        --i
                        jLayoutList.removeAt(j)
                        break
                    }
                    ++j
                }
                ++i
            }
        }

        // sort them by support min
        vdjList.sortByDescending({ vdj -> vdj.supportMin })

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

        // add all of the VDJ with only one side
        for ((geneType, layout) in oneSidedLayouts.entries())
        {
            val vdj = tryCreateOneSidedVdj(geneType, layout)
            if (vdj != null)
            {
                vdjList.add(vdj)
            }
        }

        // sort again
        vdjList.sortByDescending({ vdj -> vdj.supportMin })

        sLogger.info("{} vdj sequences after merge", vdjList.size)

        return vdjList
    }

    // This function tries to complete the layout by searching for the missing anchor
    // we already have one side of the anchor
    // use the blosum searcher to find the other side
    // if it is found then build a VDJSequence
    fun tryCompleteLayoutWithBlosum(layoutGeneType: VJGeneType, layout: ReadLayout, seqIndex: Int)
            : VDJSequence?
    {
        sLogger.debug("try complete {} layout: {}", layoutGeneType, layout.consensusSequence())

        val targetAnchorType = CiderUtils.getPairedVjGeneType(layoutGeneType)

        // we want to use the indices to work where things are
        val layoutSeq: String = layout.consensusSequence()

        val layoutAnchorRange: IntRange? = vjLayoutAdaptor.getAnchorRange(layoutGeneType.vj, layout)

        if (layoutAnchorRange == null)
            return null

        val searchStart: Int
        val searchEnd: Int

        if (targetAnchorType.vj == VJ.V)
        {
            searchStart = 0
            searchEnd = layoutAnchorRange.first
        }
        else
        {
            searchStart = layoutAnchorRange.last + 1
            searchEnd = layoutSeq.length
        }

        // find all the homolog sequences
        val anchorBlosumMatch: AnchorBlosumMatch? =
            anchorBlosumSearcher.searchForAnchor(layoutSeq, targetAnchorType, searchStart, searchEnd)

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
                geneType = targetAnchorType,
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
                geneType = targetAnchorType,
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

        sLogger.debug("built VDJ sequence: {}-{}-{}",
            Codons.aminoAcidFromBases(vdj.vAnchorSequence),
            Codons.aminoAcidFromBases(vdj.cdr3SequenceShort),
            Codons.aminoAcidFromBases(vdj.jAnchorSequence))

        return vdj
    }

    // we should probably say allow some more leeway if a position has very little support
    fun tryOverlapVJ(vLayout: ReadLayout, jLayout: ReadLayout, vLayoutGeneType: VJGeneType, jLayoutGeneType: VJGeneType): VDJSequence?
    {
        assert(vLayoutGeneType.vj == VJ.V)
        assert(jLayoutGeneType.vj == VJ.J)

        if (vLayoutGeneType != CiderUtils.getPairedVjGeneType(jLayoutGeneType))
        {
            // mismatched type, i.e. IGHV must pair with IGHJ
            return null
        }

        val vLayoutSeq = vLayout.consensusSequence()
        val jLayoutSeq = jLayout.consensusSequence()

        var overlapBases: Int = -1
        var matchedBases: Int = 0

        // try multiple location to align them
        for (i in minOverlappedBases until Math.min(vLayoutSeq.length, jLayoutSeq.length))
        {
            var seqMatch = true
            matchedBases = 0

            // check for overlap
            for (j in 0 until i)
            {
                val vIndex = vLayoutSeq.length - i + j
                val jIndex = j

                val vBase = vLayoutSeq[vIndex]
                val jBase = jLayoutSeq[jIndex]

                if (vBase == jBase)
                {
                    ++matchedBases
                }
                else if (vBase != 'N' && jBase != 'N')
                {
                    seqMatch = false
                    break
                }
            }

            if (seqMatch)
            {
                // match found! now put them together into a VDJ sequence
                overlapBases = i
                break
            }
        }

        if (overlapBases == -1)
            return null

        // put the layouts together into one
        // consider the following, V/J is the layout aligned positions, + are the overlap region:
        // -----V------+++++++++------J------
        //      |_____________________|   <---- this is the amount we need to shift the aligned position of the J layout reads
        // And this amount is vLayoutLength - V + J - overlapBases
        val alignedPositionShift = -(vLayout.length - vLayout.alignedPosition + jLayout.alignedPosition - overlapBases)

        val combinedVjLayout: ReadLayout = ReadLayout.merge(vLayout, jLayout, 0, alignedPositionShift, minBaseQuality)
        combinedVjLayout.id = "${vLayout.id},${jLayout.id}"

        var vAnchorRange: IntRange? = vjLayoutAdaptor.getAnchorRange(VJ.V, vLayout)
        var jAnchorRange: IntRange? = vjLayoutAdaptor.getAnchorRange(VJ.J, jLayout)

        if (vAnchorRange == null || jAnchorRange == null)
            return null

        sLogger.debug("overlap bases: {}, matched bases: {}, vLayout: {}, jLayout: {}, combinedLayout: {}",
            overlapBases, matchedBases, vLayout.consensusSequence(), jLayout.consensusSequence(), combinedVjLayout.consensusSequence())

        val layoutSliceStart = vAnchorRange.first

        // make sure we shift the anchor ranges as well now we trimmed the left
        vAnchorRange = vAnchorRange.first - layoutSliceStart .. vAnchorRange.last - layoutSliceStart

        // we need to shift the jAnchor range to be the range within layout
        // ---VVVVVV------+++++++++------JJJJJJ---
        //    |___________|  <---- this is the amount we need to add to the j anchor range
        // which is vLayout length - layoutSliceStart - overlap bases
        val jAnchorRangeShift = vLayout.length - layoutSliceStart - overlapBases
        jAnchorRange = jAnchorRange.first + jAnchorRangeShift .. jAnchorRange.last + jAnchorRangeShift

        val layoutSliceEnd = layoutSliceStart + jAnchorRange.last + 1 // inclusive to exclusive

        val vAnchorBoundary = vAnchorRange.last + 1
        val jAnchorBoundary = jAnchorRange.first

        if (jAnchorBoundary - vAnchorBoundary < CiderConstants.MIN_CDR3_LENGTH_BASES)
            return null

        // construct a VDJ sequence
        val vAnchor: VJAnchor
        val jAnchor: VJAnchor

        vAnchor = createVJAnchorByReadMatch(
            anchorBoundary = vAnchorBoundary,
            vjGeneType = vLayoutGeneType,
            layout = vLayout)
        jAnchor = createVJAnchorByReadMatch(
            anchorBoundary = jAnchorBoundary,
            vjGeneType = jLayoutGeneType,
            layout = jLayout)

        val vdj = VDJSequence(combinedVjLayout, layoutSliceStart, layoutSliceEnd, vAnchor, jAnchor)

        sLogger.debug("built VDJ sequence: {}-{}-{}, by overlapping V layout({}): {} and J layout({}): {}",
            Codons.aminoAcidFromBases(vdj.vAnchorSequence),
            Codons.aminoAcidFromBases(vdj.cdr3SequenceShort),
            Codons.aminoAcidFromBases(vdj.jAnchorSequence),
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
            val layoutEnd = layout.length

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

            val layoutStart = 0
            val layoutEnd = jAnchorRange.last + 1

            jAnchor = createVJAnchorByReadMatch(
                anchorBoundary = jAnchorRange.first,
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

                if (vdjSequenceCompare(vdj1, vdj2, diffAccumulator))
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
            require(vdj1.vAnchor != null)
            require(vdj1.jAnchor != null)
            require(vdj2.vAnchor != null)
            require(vdj2.jAnchor != null)

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
                vdj1.layoutSliceStart, vdj1.layoutSliceEnd, vdj1.vAnchor.anchorBoundary, vdj1.jAnchor.anchorBoundary,
                vdj2.aminoAcidSequenceFormatted, vdj2.vAnchor.matchMethod, vdj2.jAnchor.matchMethod, vdj2.layout.alignedPosition,
                vdj2.layoutSliceStart, vdj2.layoutSliceEnd, vdj2.vAnchor.anchorBoundary, vdj2.jAnchor.anchorBoundary)

            val vAnchorBoundary = Math.max(vdj1.vAnchor.anchorBoundary, vdj2.vAnchor.anchorBoundary)
            val jAnchorBoundary = Math.max(vdj1.jAnchor.anchorBoundary, vdj2.jAnchor.anchorBoundary)

            // merge the anchor objects together
            val vAnchor: VJAnchor = mergeAnchors(vdj1.vAnchor, vdj2.vAnchor, vdj1.vAnchorLength, vdj2.vAnchorLength, vAnchorBoundary)
            val jAnchor: VJAnchor = mergeAnchors(vdj1.jAnchor, vdj2.jAnchor, vdj1.jAnchorLength, vdj2.jAnchorLength, jAnchorBoundary)

            val jAnchorLength: Int = Math.max(vdj1.jAnchorLength, vdj2.jAnchorLength)

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

            require(vdjPrimary.vAnchor != null)
            require(vdjPrimary.jAnchor != null)
            require(vdjSecondary.vAnchor != null)
            require(vdjSecondary.jAnchor != null)

            // we calculate the aligned position for both sequences with respect to the V anchor boundary
            val primaryAlignedPos = vdjPrimary.layout.alignedPosition - vdjPrimary.layoutSliceStart - vdjPrimary.vAnchor.anchorBoundary
            val secondaryAlignedPos = vdjSecondary.layout.alignedPosition - vdjSecondary.layoutSliceStart - vdjSecondary.vAnchor.anchorBoundary

            // this would shift them to overlap each other
            val secondaryLayoutOffsetShift = primaryAlignedPos - secondaryAlignedPos
            val mergedLayout = ReadLayout.merge(vdjPrimary.layout, vdjSecondary.layout, 0, secondaryLayoutOffsetShift, minBaseQuality)
            mergedLayout.id = "${vdjPrimary.layout.id};${vdjSecondary.layout.id}"

            // now we have to work out where to place the layout, using the final aligned position
            val posStartWithinLayout = vdjPrimary.layoutSliceStart +
                    vdjPrimary.vAnchor.anchorBoundary - vAnchor.anchorBoundary + // adjust for V anchor length
                    mergedLayout.alignedPosition - vdjPrimary.layout.alignedPosition // adjust for align position difference
            val posEndWithinLayout = posStartWithinLayout + jAnchor.anchorBoundary + jAnchorLength

            sLogger.debug("primary start: {}, primary v anchor: {}, v anchor: {}, layout aligned pos: {}, primary aligned pos: {}",
                vdjPrimary.layoutSliceStart, vdjPrimary.vAnchor.anchorBoundary, vAnchor.anchorBoundary,
                mergedLayout.alignedPosition, vdjPrimary.layout.alignedPosition)

            sLogger.debug("layout: {}, aligned pos: {}, pos within layout: {}-{}, v: {}, j: {}",
                mergedLayout.consensusSequence(), mergedLayout.alignedPosition, posStartWithinLayout, posEndWithinLayout,
                vAnchor.anchorBoundary, jAnchor.anchorBoundary)

            /*assert(supports.size > jAnchor.anchorBoundary)

            if (supports.size < jAnchorBoundary)
                throw RuntimeException()*/

            val combinedVDJ = VDJSequence(mergedLayout, posStartWithinLayout, posEndWithinLayout, vAnchor, jAnchor)

            require(combinedVDJ.vAnchor != null)
            require(combinedVDJ.jAnchor != null)

            sLogger.debug("merge: {}(v:{}, j:{}, read count: {}) and {}(v:{}, j:{}, read count: {}) together, merged: {}(v:{}, j:{}), 2ndary aligned pos shift: {}",
                vdj1.aminoAcidSequenceFormatted, vdj1.vAnchor.matchMethod, vdj1.jAnchor.matchMethod, vdj1.numReads,
                vdj2.aminoAcidSequenceFormatted, vdj2.vAnchor.matchMethod, vdj2.jAnchor.matchMethod, vdj2.numReads,
                combinedVDJ.aminoAcidSequenceFormatted, combinedVDJ.vAnchor.matchMethod, combinedVDJ.jAnchor.matchMethod,
                secondaryLayoutOffsetShift)

            return combinedVDJ
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
            return vdjSequenceCompare(vdj1, vdj2, diffAccumulator)
        }

        // TODO write unit test
        private fun vdjSequenceCompare(vdj1: VDJSequence, vdj2: VDJSequence,
                                       diffAccumulator: (Map.Entry<Char, Int>, Map.Entry<Char, Int>) -> Boolean) : Boolean
        {
            require(vdj1.vAnchor != null)
            require(vdj1.jAnchor != null)
            require(vdj2.vAnchor != null)
            require(vdj2.jAnchor != null)

            // the lengths between V and J anchor must be equal
            if (vdj1.jAnchor.anchorBoundary - vdj1.vAnchor.anchorBoundary !=
                vdj2.jAnchor.anchorBoundary - vdj2.vAnchor.anchorBoundary)
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
            require(vdj1.vAnchor != null)
            require(vdj1.jAnchor != null)
            require(vdj2.vAnchor != null)
            require(vdj2.jAnchor != null)

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

        fun mergeAnchors(anchor1: VJAnchor, anchor2: VJAnchor,
                         anchor1Length: Int, anchor2Length: Int,
                         newAnchorBoundary: Int) : VJAnchor
        {
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
            val anchorTemplate: VJAnchor = if (anchor1Length >= anchor2Length) anchor1 else anchor2

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

            throw IllegalArgumentException("unknown VJAnchor subtype: ${anchorTemplate::class}")
        }
    }
}