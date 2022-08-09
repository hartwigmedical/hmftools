package com.hartwig.hmftools.cdr3

import com.google.common.collect.ArrayListMultimap
import com.hartwig.hmftools.cdr3.layout.ReadLayout
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
        val orphanedLayouts: ArrayListMultimap<VJGeneType, ReadLayout> = ArrayListMultimap.create()

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
                    orphanedLayouts.put(geneType, l)
                }

                ++seqIndex
            }
        }

        // try to join together the orphaned layouts
        for (vGeneType: VJGeneType in orphanedLayouts.keySet().filter({ vjGeneType -> vjGeneType.vj == VJ.V }))
        {
            val jGeneType: VJGeneType = Cdr3Utils.getPairedVjGeneType(vGeneType)
            val vLayoutList: List<ReadLayout> = orphanedLayouts[vGeneType]
            val jLayoutList: MutableList<ReadLayout> = orphanedLayouts[jGeneType]

            // we try each pair see if they can be aligned together
            // we find the best aligned pair and remove them and try again
            for (vLayout in vLayoutList)
            {
                var j = 0
                while (j < jLayoutList.size)
                {
                    val jLayout: ReadLayout = jLayoutList[j]
                    val vdj: VDJSequence? = tryOverlapVJ(vLayout, jLayout, vGeneType, jGeneType)

                    if (vdj != null)
                    {
                        vdjList.add(vdj)
                        jLayoutList.removeAt(j)
                        break
                    }
                    ++j
                }
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

    // This function tries to complete the layout by searching for the missing anchor
    // we already have one side of the anchor
    // use the blosum searcher to find the other side
    // if it is found then build a VDJSequence
    fun tryCompleteLayoutWithBlosum(layoutGeneType: VJGeneType, layout: ReadLayout, seqIndex: Int)
        : VDJSequence?
    {
        val targetAnchorType = Cdr3Utils.getPairedVjGeneType(layoutGeneType)

        // we want to use the indices to work where things are
        val layoutSeq: String = layout.consensusSequence()

        val layoutAnchorRange: IntRange = vjLayoutAdaptor.getAnchorRange(layoutGeneType.vj, layout)!!

        val searchStart = if (targetAnchorType.vj == VJ.V) 0 else layoutAnchorRange.last + 1
        val searchEnd = if (targetAnchorType.vj == VJ.V) layoutAnchorRange.first else layoutSeq.length

        // find all the homolog sequences
        val anchorBlosumMatch: AnchorBlosumMatch? =
            anchorBlosumSearcher.searchForAnchor(layoutSeq, targetAnchorType, searchStart, searchEnd)

        if (anchorBlosumMatch == null)
            return null

        var vAnchorRange: IntRange
        var jAnchorRange: IntRange

        if (layoutGeneType.vj == VJ.V)
        {
            vAnchorRange = vjLayoutAdaptor.getAnchorRange(layoutGeneType.vj, layout)!!
            jAnchorRange = anchorBlosumMatch.anchorStart until anchorBlosumMatch.anchorEnd
        }
        else
        {
            vAnchorRange = anchorBlosumMatch.anchorStart until anchorBlosumMatch.anchorEnd
            jAnchorRange = vjLayoutAdaptor.getAnchorRange(layoutGeneType.vj, layout)!!
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
                vjGeneType = layoutGeneType,
                matchMethod = vjLayoutAdaptor.getAnchorMatchType(layout)
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
                vjGeneType = layoutGeneType,
                matchMethod = vjLayoutAdaptor.getAnchorMatchType(layout)
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
    fun tryOverlapVJ(vLayout: ReadLayout, jLayout: ReadLayout, vLayoutGeneType: VJGeneType, jLayoutGeneType: VJGeneType): VDJSequence?
    {
        assert(vLayoutGeneType.vj == VJ.V)
        assert(jLayoutGeneType.vj == VJ.J)

        if (vLayoutGeneType != Cdr3Utils.getPairedVjGeneType(jLayoutGeneType))
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

        // construct a VDJ sequence
        val vAnchor: VJAnchor
        val jAnchor: VJAnchor

        vAnchor = createVJAnchorByReadMatch(
            type = VJAnchor.Type.V,
            anchorBoundary = vAnchorRange.last + 1,
            vjGeneType = vLayoutGeneType,
            matchMethod = vjLayoutAdaptor.getAnchorMatchType(vLayout))
        jAnchor = createVJAnchorByReadMatch(
            type = VJAnchor.Type.J,
            anchorBoundary = jAnchorRange.first,
            vjGeneType = jLayoutGeneType,
            matchMethod = vjLayoutAdaptor.getAnchorMatchType(jLayout))

        val vdj = VDJSequence("", combinedVjLayout, layoutSliceStart, layoutSliceEnd, vAnchor, jAnchor)

        sLogger.debug("built VDJ sequence: {}-{}-{}, by overlapping V layout({}): {} and J layout({}): {}",
            Codons.aminoAcidFromBases(vdj.vAnchorSequence),
            Codons.aminoAcidFromBases(vdj.dSequenceShort),
            Codons.aminoAcidFromBases(vdj.jAnchorSequence),
            vLayout.id, vLayout.consensusSequence(), jLayout.id, jLayout.consensusSequence())

        return vdj
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
            val primaryAlignedPos = vdjPrimary.layout.alignedPosition - vdjPrimary.layoutSliceStart - vdjPrimary.vAnchor.anchorBoundary
            val secondaryAlignedPos = vdjSecondary.layout.alignedPosition - vdjSecondary.layoutSliceStart - vdjSecondary.vAnchor.anchorBoundary

            // this would shift them to overlap each other
            val secondaryLayoutOffsetShift = primaryAlignedPos - secondaryAlignedPos
            val mergedLayout = ReadLayout.merge(vdjPrimary.layout, vdjSecondary.layout, 0, secondaryLayoutOffsetShift, minBaseQuality)
            mergedLayout.id = "${vdjPrimary.id},${vdjSecondary.id}"

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

            val combinedVDJ = VDJSequence(vdjPrimary.id, mergedLayout, posStartWithinLayout, posEndWithinLayout, vAnchor, jAnchor)

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

        fun copyVjAnchorTemplate(anchorTemplate: VJAnchor, anchorBoundary: Int): VJAnchor
        {
            if (anchorTemplate is VJAnchorByReadMatch)
            {
                // smart casted to VJAnchorByReadMatch
                return anchorTemplate.copy(anchorBoundary = anchorBoundary)
            }
            else if (anchorTemplate is VJAnchorByBlosum)
            {
                // smart casted to VJAnchorByBlosum
                return anchorTemplate.copy(anchorBoundary = anchorBoundary)
            }
            else
            {
                throw IllegalArgumentException("unknown VJAnchor child type")
            }
        }

        fun createVJAnchorByReadMatch(
            type: VJAnchor.Type,
            anchorBoundary: Int,
            vjGeneType: VJGeneType,
            matchMethod: VJReadCandidate.AnchorMatchMethod) : VJAnchorByReadMatch
        {
            val matchMethodStr =
                when (matchMethod)
                {
                    VJReadCandidate.AnchorMatchMethod.ALIGN -> "align"
                    VJReadCandidate.AnchorMatchMethod.EXACT -> "exact"
                }
            return VJAnchorByReadMatch(type, vjGeneType, anchorBoundary, matchMethodStr)
        }
    }
}