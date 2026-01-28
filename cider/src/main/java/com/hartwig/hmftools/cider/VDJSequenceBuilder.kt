package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.CiderConstants.LAYOUT_MIN_READ_OVERLAP_BASES
import com.hartwig.hmftools.cider.CiderConstants.LAYOUT_MATE_EXTEND_MISMATCHES_MAX
import com.hartwig.hmftools.cider.CiderConstants.MIN_ANCHOR_LENGTH_BASES
import com.hartwig.hmftools.cider.genes.VJ
import com.hartwig.hmftools.cider.genes.VJGeneType
import com.hartwig.hmftools.cider.layout.ReadLayout
import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.perf.TaskExecutor
import htsjdk.samtools.util.SequenceUtil.N
import org.apache.logging.log4j.LogManager
import java.util.*
import kotlin.collections.ArrayList
import kotlin.math.max
import kotlin.math.min

class VDJSequenceBuilder(private val vjLayoutAdaptor: IVJReadLayoutAdaptor,
                         private val anchorBlosumSearcher: IAnchorBlosumSearcher,
                         private val minBaseQuality: Byte,
                         private val minOverlappedBases: Int)
{
    fun buildVDJSequences(
        layoutMultimap: Map<VJGeneType, List<ReadLayout>>, mateLayoutReads: Map<String, List<ReadLayout.Read>>, threadCount: Int
    ) : List<VDJSequence>
    {
        // Care is taken here to populate this list in a deterministic order, otherwise the output can be somewhat different between runs,
        // which is confusing.
        // Note that the iteration order of maps and sets is not defined and may vary.
        val vdjList = ArrayList<VDJSequence>()

        // remember the layouts which are not used
        val oneSidedLayouts = HashMap<VJGeneType, ArrayList<ReadLayout>>()
        for ((geneType, layouts) in layoutMultimap.entries)
        {
            oneSidedLayouts[geneType] = ArrayList(layouts)
        }

        // try to join together the layouts
        joinVjLayouts(oneSidedLayouts, mateLayoutReads, vdjList)
        sLogger.info("built {} vdj by joining V with J layouts", vdjList.size)

        // for the remaining ones, we try to complete with blosum
        for ((geneType, layouts) in oneSidedLayouts.entries.sortedBy { it.key })
        {
            val itr = layouts.iterator()
            while (itr.hasNext())
            {
                val layout = itr.next()
                val vdj = tryCompleteLayoutWithBlosum(geneType, layout)

                if (vdj != null)
                {
                    vdjList.add(vdj)
                    itr.remove()
                }
            }
        }

        // add all of the VDJ with only one side
        for ((geneType, layouts) in oneSidedLayouts.entries.sortedBy { it.key })
        {
            for (layout in layouts)
            {
                val vdj = tryCreateOneSidedVdj(geneType, layout)
                if (vdj != null)
                {
                    vdjList.add(vdj)
                }
            }
        }

        // sort them by number of reads
        vdjList.sortByDescending { vdj -> vdj.numReads }

        sLogger.info("found {} vdj sequences before merge", vdjList.size)

        val vjGeneToVdjMap = HashMap<Pair<VJGeneType?, VJGeneType?>, ArrayList<VDJSequence>>()
        for (vdj in vdjList)
        {
            if (vdj.vAnchor == null && vdj.jAnchor == null)
            {
                assert(false)
                continue
            }
            val key = Pair(vdj.vAnchor?.geneType, vdj.jAnchor?.geneType)
            vjGeneToVdjMap.computeIfAbsent(key) { ArrayList() }.add(vdj)
        }

        vdjList.clear()
        vdjList.addAll(mergeIdentical(vjGeneToVdjMap, minBaseQuality, threadCount))

        // sort again
        vdjList.sortByDescending { vdj -> vdj.numReads }

        sLogger.info("{} vdj sequences after merge", vdjList.size)

        return vdjList
    }

    private fun joinVjLayouts(
        oneSidedLayouts: Map<VJGeneType, MutableList<ReadLayout>>, mateLayoutReads: Map<String, List<ReadLayout.Read>>,
        vdjList: MutableList<VDJSequence>)
    {
        val geneTypes = oneSidedLayouts.keys.filter { vjGeneType -> vjGeneType.vj == VJ.V }.sorted()
        for (vGeneType: VJGeneType in geneTypes)
        {
            for (jGeneType: VJGeneType in vGeneType.pairedVjGeneTypes())
            {
                joinVjLayouts(vGeneType, jGeneType,
                    oneSidedLayouts[vGeneType] ?: ArrayList(),
                    oneSidedLayouts[jGeneType] ?: ArrayList(),
                    mateLayoutReads,
                    vdjList)
            }
        }
    }

    private fun joinVjLayouts(
        vGeneType: VJGeneType, jGeneType: VJGeneType, vLayouts: MutableList<ReadLayout>, jLayouts: MutableList<ReadLayout>,
        mateLayoutReads: Map<String, List<ReadLayout.Read>>, vdjList: MutableList<VDJSequence>)
    {
        joinVjLayoutsBySharedReads(vGeneType, jGeneType, vLayouts, jLayouts, mateLayoutReads, vdjList)
        joinVjLayoutsByWordHash(vGeneType, jGeneType, vLayouts, jLayouts, mateLayoutReads, vdjList)
    }

    // this method of joining V and J layout finds the candidate matches by looking for layouts
    // that share reads. The layout that share the largest number of reads are tested first.
    private fun joinVjLayoutsBySharedReads(
        vGeneType: VJGeneType, jGeneType: VJGeneType, vLayouts: MutableList<ReadLayout>, jLayouts: MutableList<ReadLayout>,
        mateLayoutReads: Map<String, List<ReadLayout.Read>>, vdjList: MutableList<VDJSequence>)
    {
        val consumedLayouts: MutableSet<ReadLayout> = Collections.newSetFromMap(IdentityHashMap())

        // first we use the read IDs
        val readIdMultimap = HashMap<String, ArrayList<ReadLayout>>()

        // put all v layout read IDs into the map
        for (jLayout in jLayouts)
        {
            jLayout.reads.forEach { r -> readIdMultimap.computeIfAbsent(r.readKey.readName) { ArrayList() }.add(jLayout) }
        }

        val vLayoutItr: MutableIterator<ReadLayout> = vLayouts.iterator()
        while (vLayoutItr.hasNext())
        {
            val vLayout: ReadLayout = vLayoutItr.next()

            // store a hash table with the layout and the number of matching reads
            // the highest matching read is tried first
            val jReadMatchCount = IdentityHashMap<ReadLayout, Int>()

            // look up using the read ids
            for (read in vLayout.reads)
            {
                (readIdMultimap[read.readKey.readName] ?: continue).forEach { l ->
                    jReadMatchCount.merge(l, 1) { exitingCount: Int, _ -> exitingCount + 1 }
                }
            }

            // sort the layouts by their number of reads that are shared
            val jCandidates = jReadMatchCount.entries
                .filter { (jLayout, _) -> !consumedLayouts.contains(jLayout) }
                // Include stable ID to produce deterministic order in the case of tied number of reads.
                .sortedWith(compareBy({ -it.value }, { it.key.id }))

            for ((jLayout, _) in jCandidates)
            {
                val vdj: VDJSequence? = tryOverlapVJ(vLayout, jLayout, vGeneType, jGeneType, mateLayoutReads)

                if (vdj != null)
                {
                    sLogger.trace("vj joined by matching reads, num candidates={}", jCandidates.size)
                    vdjList.add(vdj)
                    vLayoutItr.remove()
                    consumedLayouts.add(jLayout)
                    break
                }
            }
        }

        // now remove all consumed layouts from the lists
        val jLayoutItr = jLayouts.iterator()
        while (jLayoutItr.hasNext())
        {
            if (consumedLayouts.contains(jLayoutItr.next()))
            {
                jLayoutItr.remove()
            }
        }
    }

    // this method of joining V and J layouts generates word hashes of each 8-mers. For V layout the
    // hash is generated from the V anchor end onwards, and for J layout the hashes are generated
    // until the J anchor start. The reason for doing it this way is that the V and J anchor sequences
    // are very repetitive and cause too many spurious hash matches (80% of j layouts will match with any
    // v layouts from my tests)
    internal fun joinVjLayoutsByWordHash(
        vGeneType: VJGeneType, jGeneType: VJGeneType, vLayouts: MutableList<ReadLayout>, jLayouts: MutableList<ReadLayout>,
        mateLayoutReads: Map<String, List<ReadLayout.Read>>, vdjList: MutableList<VDJSequence>)
    {
        // for v and j layouts that have not found a match yet,
        // we calculate a key for each 5 high qual bases. The hash is actually a match since
        // we use 2bits per base
        val baseHashMultimap = HashMap<Int, ArrayList<ReadLayout>>()

        for (jLayout in jLayouts)
        {
            val layoutAnchorRange: IntRange? = vjLayoutAdaptor.getAnchorRange(VJ.J, jLayout)

            if (layoutAnchorRange == null || layoutAnchorRange.last < 0)
                continue

            var seq = String(jLayout.highConfidenceSequence(CiderConstants.MIN_VJ_LAYOUT_HIGH_QUAL_READ_FRACTION))

            // for j layout, we only hash bases up to the anchor start, reason is that the anchor sequences
            // are too repetitive to make hashing meaningful
            seq = seq.substring(0, min(layoutAnchorRange.first, seq.length))

            val hashList = VdjBuilderUtils.calcSequenceWordHashes(seq, CiderConstants.VJ_JOIN_HASH_WORD_SIZE)
            for (hash in hashList)
            {
                baseHashMultimap.computeIfAbsent(hash) { ArrayList() }.add(jLayout)
            }
        }

        val consumedLayouts: MutableSet<ReadLayout> = Collections.newSetFromMap(IdentityHashMap())

        // now we got them all hashed, we can find the the v layout matches
        val vLayoutItr: MutableIterator<ReadLayout> = vLayouts.iterator()
        while (vLayoutItr.hasNext())
        {
            val vLayout: ReadLayout = vLayoutItr.next()

            val layoutAnchorRange: IntRange? = vjLayoutAdaptor.getAnchorRange(VJ.V, vLayout)

            if (layoutAnchorRange == null || layoutAnchorRange.last < 0)
                continue

            var seq = String(vLayout.highConfidenceSequence(CiderConstants.MIN_VJ_LAYOUT_HIGH_QUAL_READ_FRACTION))

            // for v layout, we only hash from anchor end, reason is that the anchor sequences
            // are too repetitive to make hashing meaningful
            seq = seq.substring(max(layoutAnchorRange.last + 1, 0), seq.length)

            val hashList = VdjBuilderUtils.calcSequenceWordHashes(seq, CiderConstants.VJ_JOIN_HASH_WORD_SIZE)

            // store a hash table with the layout and the number of matching reads
            // the highest matching read is tried first
            val jHashMatchCount = IdentityHashMap<ReadLayout, Int>()

            // look up using the read ids
            for (hash in hashList)
            {
                (baseHashMultimap[hash] ?: continue).forEach { l ->
                    jHashMatchCount.merge(l, 1) { exitingCount: Int, _ -> exitingCount + 1 }
                }
            }

            // sort descending by the number of hashes that are shared
            val jCandidates = jHashMatchCount.entries
                .filter { (jLayout, _) -> !consumedLayouts.contains(jLayout) }
                // Include stable ID to produce deterministic order in the case of tied number of hashes.
                .sortedWith(compareBy({ -it.value }, { it.key.id }))

            for ((jLayout, _) in jCandidates)
            {
                val vdj: VDJSequence? = tryOverlapVJ(vLayout, jLayout, vGeneType, jGeneType, mateLayoutReads)

                if (vdj != null)
                {
                    sLogger.trace("vj joined by hash, num candidates={}", jCandidates.size)
                    vdjList.add(vdj)
                    vLayoutItr.remove()
                    consumedLayouts.add(jLayout)
                    break
                }
            }
        }

        // now remove all consumed layouts from the j layout list
        val jLayoutItr = jLayouts.iterator()
        while (jLayoutItr.hasNext())
        {
            if (consumedLayouts.contains(jLayoutItr.next()))
            {
                jLayoutItr.remove()
            }
        }
    }

    // This function tries to complete the layout by searching for the missing anchor
    // we already have one side of the anchor
    // use the blosum searcher to find the other side
    // if it is found then build a VDJSequence
    internal fun tryCompleteLayoutWithBlosum(layoutGeneType: VJGeneType, layout: ReadLayout)
            : VDJSequence?
    {
        sLogger.trace("try complete {} layout: {}", layoutGeneType, layout.consensusSequenceString())

        val targetAnchorTypes = layoutGeneType.pairedVjGeneTypes()

        // we want to use the indices to work where things are
        val layoutSeq: String = layout.consensusSequenceString()

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

        var vAnchorRange: IntRange
        var jAnchorRange: IntRange

        if (layoutGeneType.vj == VJ.V)
        {
            vAnchorRange = layoutAnchorRange
            jAnchorRange = anchorBlosumMatch.anchorStart until anchorBlosumMatch.anchorEnd
        }
        else
        {
            vAnchorRange = anchorBlosumMatch.anchorStart until anchorBlosumMatch.anchorEnd
            jAnchorRange = layoutAnchorRange
        }

        val layoutStart = max(vAnchorRange.first, 0)
        val layoutEnd = min(jAnchorRange.last + 1, layout.length) // inclusive to exclusive

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

        sLogger.trace("built VDJ sequence: {}, {}",
            vdj.aminoAcidSequenceFormatted, vdj.sequenceFormatted)

        return vdj
    }

    // we should probably say allow some more leeway if a position has very little support
    internal fun tryOverlapVJ(
        vLayout: ReadLayout, jLayout: ReadLayout, vLayoutGeneType: VJGeneType, jLayoutGeneType: VJGeneType,
        mateLayoutReads: Map<String, List<ReadLayout.Read>>): VDJSequence?
    {
        assert(vLayoutGeneType.vj == VJ.V)
        assert(jLayoutGeneType.vj == VJ.J)

        if (!vLayoutGeneType.pairedVjGeneTypes().contains(jLayoutGeneType))
        {
            // mismatched type, i.e. IGHV must pair with IGHJ
            return null
        }

        val vLayoutSeq = vLayout.highConfidenceSequence(CiderConstants.MIN_VJ_LAYOUT_HIGH_QUAL_READ_FRACTION)
        val jLayoutSeq = jLayout.highConfidenceSequence(CiderConstants.MIN_VJ_LAYOUT_HIGH_QUAL_READ_FRACTION)

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

        val combinedVjLayout: ReadLayout = ReadLayout.merge(vLayout, jLayout, 0, jAlignedPositionShift, minBaseQuality)
        combinedVjLayout.id = "${vLayout.id};${jLayout.id}"

        // Extend the layout sequence by adding mates of reads in the layout.
        // These mates wouldn't have been included earlier because they are not close enough to the anchor.
        tryExtendCombinedLayoutWithMates(combinedVjLayout, mateLayoutReads, vLayoutGeneType, jLayoutGeneType, jAlignedPositionShift)

        // next we want to work out where the anchor ranges are in the combined layout
        // we use "extrapolated" anchor range, the reason is that the layout might contain only half
        // of the anchor, if we do not use extrapolated anchor range then the final anchor length
        // could be wrong
        // for example:
        // v layout:   TGC-GAATACC-CACATCCTGA-GAGTGG-TCAGATA
        // j layout:                              GG-TCAGATAACT
        //                 |_____|            |____|
        //            V anchor(align)    J anchor(align)
        // here the J layout only contains the first 2 bases of the J anchor. Using extrapolated anchor
        // range, we get anchor range of (-4, 2), and we can get the full j anchor in the merged
        // layout.
        var vAnchorRange: IntRange = vjLayoutAdaptor.getExtrapolatedAnchorRange(VJ.V, vLayout)
        var jAnchorRange: IntRange = vjLayoutAdaptor.getExtrapolatedAnchorRange(VJ.J, jLayout)

        // now we use the aligned position of the combined layout
        val vAnchorShift = combinedVjLayout.alignedPosition - vLayout.alignedPosition
        vAnchorRange = vAnchorRange.first + vAnchorShift .. vAnchorRange.last + vAnchorShift

        val jAnchorShift = combinedVjLayout.alignedPosition - jLayout.alignedPosition - jAlignedPositionShift
        jAnchorRange = jAnchorRange.first + jAnchorShift .. jAnchorRange.last + jAnchorShift

        sLogger.trace("overlap: {}, jAlignShift: {}, vAnchor: {}, jAnchor: {}, vLayout: {}, jLayout: {}, combinedLayout: {}",
            vjLayoutOverlap, jAlignedPositionShift, vAnchorRange, jAnchorRange,
            vLayout.consensusSequenceString(), jLayout.consensusSequenceString(), combinedVjLayout.consensusSequenceString())

        val layoutSliceStart: Int
        val layoutSliceEnd = min(jAnchorRange.last + 1, combinedVjLayout.length) // inclusive to exclusive

        var vAnchorBoundary = vAnchorRange.last + 1
        var jAnchorBoundary = jAnchorRange.first

        // we need to have at least 3 bases in the anchor
        if (vAnchorBoundary < MIN_ANCHOR_LENGTH_BASES || vAnchorBoundary > combinedVjLayout.length - MIN_ANCHOR_LENGTH_BASES ||
            jAnchorBoundary < MIN_ANCHOR_LENGTH_BASES || jAnchorBoundary > combinedVjLayout.length - MIN_ANCHOR_LENGTH_BASES)
            return null

        // construct a VDJ sequence
        val vAnchor: VJAnchor?
        val jAnchor: VJAnchor

        //if (jAnchorBoundary - vAnchorBoundary >= 0)
        if (vAnchorRange.last < jAnchorRange.first)
        {
            layoutSliceStart = max(vAnchorRange.first, 0)
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
            sLogger.trace("VDJ vAnchorBoundary({}) < jAnchorBoundary({}), setting v anchor to null",
                vAnchorBoundary, jAnchorBoundary)
        }

        jAnchor = createVJAnchorByReadMatch(
            anchorBoundary = jAnchorBoundary,
            vjGeneType = jLayoutGeneType,
            layout = jLayout)

        val vdj = VDJSequence(combinedVjLayout, layoutSliceStart, layoutSliceEnd, vAnchor, jAnchor)

        sLogger.trace("built VDJ sequence: {}, by overlapping V layout({}): {} and J layout({}): {}",
            vdj.aminoAcidSequenceFormatted,
            vLayout.id, vLayout.consensusSequenceString(), jLayout.id, jLayout.consensusSequenceString())

        return vdj
    }

    internal fun tryCreateOneSidedVdj(layoutGeneType: VJGeneType, layout: ReadLayout): VDJSequence?
    {
        sLogger.trace("create one sided {} layout: {}", layoutGeneType, layout.consensusSequenceString())

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
            val vAnchorRange = layoutAnchorRange

            val layoutStart = max(vAnchorRange.first, 0)

            // we limit one sided VDJ post V anchor length
            val layoutEnd = min(layout.length, vAnchorRange.last + 1 + CiderConstants.PARTIAL_VDJ_UNANCHORED_LENGTH_BASES)

            vAnchor = createVJAnchorByReadMatch(
                anchorBoundary = vAnchorRange.last + 1 - layoutStart,
                vjGeneType = layoutGeneType,
                layout = layout
            )
            vdj = VDJSequence(layout, layoutStart, layoutEnd, vAnchor, null)

            sLogger.trace("built V only sequence: {}-{}",
                Codons.aminoAcidFromBases(vdj.vAnchorSequence),
                Codons.aminoAcidFromBases(vdj.cdr3SequenceShort))
        }
        else
        {
            val jAnchorRange = layoutAnchorRange

            // we limit one sided VDJ pre J anchor length
            val layoutStart = max(0, jAnchorRange.first - CiderConstants.PARTIAL_VDJ_UNANCHORED_LENGTH_BASES)
            val layoutEnd = min(layout.length, jAnchorRange.last + 1)

            jAnchor = createVJAnchorByReadMatch(
                anchorBoundary = jAnchorRange.first - layoutStart,
                vjGeneType = layoutGeneType,
                layout = layout
            )
            vdj = VDJSequence(layout, layoutStart, layoutEnd, null, jAnchor)

            sLogger.trace("built J only sequence: {}-{}",
                Codons.aminoAcidFromBases(vdj.cdr3SequenceShort),
                Codons.aminoAcidFromBases(vdj.jAnchorSequence))
        }

        return vdj
    }

    private fun createVJAnchorByReadMatch(
        anchorBoundary: Int,
        layout: ReadLayout,
        vjGeneType: VJGeneType) : VJAnchorByReadMatch
    {
        val matchMethodStr =
            when (vjLayoutAdaptor.getAnchorMatchMethod(layout))
            {
                VJReadCandidate.MatchMethod.ALIGN -> "align"
                VJReadCandidate.MatchMethod.BLOSUM -> "blosum"
            }

        val templateAnchorSeq = vjLayoutAdaptor.getTemplateAnchorSequence(layout)

        return VJAnchorByReadMatch(vjGeneType.vj, vjGeneType, anchorBoundary, matchMethodStr,
            templateAnchorSeq, layout.reads.size)
    }

    private fun tryExtendCombinedLayoutWithMates(
        combinedVjLayout: ReadLayout, mateLayoutReads: Map<String, List<ReadLayout.Read>>,
        vLayoutGeneType: VJGeneType, jLayoutGeneType: VJGeneType, jAlignedPositionShift: Int)
    {
        val mates = combinedVjLayout.reads
            .flatMap { mateLayoutReads[it.readKey.readName] ?: emptyList() }
            .mapNotNull {
                val mateCandidate = vjLayoutAdaptor.toReadCandidate(it)
                val canExtend = mateCandidate.vjGeneType == vLayoutGeneType || mateCandidate.vjGeneType == jLayoutGeneType
                if (canExtend)
                {
                    // The combined layout is relative to the V side. So if the read is for the J side, need to adjust its aligned position.
                    if (mateCandidate.vjGeneType.vj == VJ.J) it.copy(it.alignedPosition + jAlignedPositionShift) else it
                }
                else null
            }
        tryExtendLayoutBySequence(combinedVjLayout, mates)
    }

    private fun tryExtendLayoutBySequence(layout: ReadLayout, reads: List<ReadLayout.Read>)
    {
        // Iteratively: collect reads which overlap the layout (by aligned position), check their sequence, and add them to the layout if
        // they mostly match.

        val remainingReads = reads.indices.toMutableSet()
        while (true)
        {
            //  layout aligned pos:              |
            //  consensus:        ----------------------------------
            //  reads: V       ----------------  |
            //         V ------------------      |
            //         J                         |          -----------------
            //         J                         |               ----------------
            //         J                         |                            -------------
            // Aligned position is the last base of the V anchor.

            // Consider all reads which overlap the current assembled sequence.
            val candidateReads = remainingReads.mapNotNull { readIndex ->
                val read = reads[readIndex]
                val readStartInLayout = layout.alignedPosition - read.alignedPosition
                val readEndInLayout = readStartInLayout + read.readLength
                val readLayoutOverlap = min(readEndInLayout, layout.length) - max(readStartInLayout, 0)
                if (readLayoutOverlap >= LAYOUT_MIN_READ_OVERLAP_BASES) Pair(readIndex, readStartInLayout) else null
            }

            // Iterate over the sequences of the candidate reads to determine if they sufficiently match the assembly.
            val candidateReadMismatches = Array(candidateReads.size) { 0 }
            for ((position, consensusBase) in layout.consensusSequence().withIndex())
            {
                for (candidateReadIndex in candidateReads.indices)
                {
                    val (readIndex, readStartInLayout) = candidateReads[candidateReadIndex]
                    val readSequence = reads[readIndex].sequence
                    val readPosition = position - readStartInLayout
                    if (readPosition >= 0 && readPosition < readSequence.size)
                    {
                        val readBase = readSequence[readPosition]
                        val basesMatch = readBase == consensusBase || readBase == N
                        if (!basesMatch)
                        {
                            candidateReadMismatches[candidateReadIndex]++
                        }
                    }
                }
            }

            val matchingReadIndices = candidateReads.indices
                .filter { candidateReadMismatches[it] <= LAYOUT_MATE_EXTEND_MISMATCHES_MAX }
                .map { candidateReads[it].first }
            if (matchingReadIndices.isEmpty())
            {
                // No more reads can extend the sequence.
                break
            }
            else
            {
                for (readIndex in matchingReadIndices)
                {
                    remainingReads.remove(readIndex)
                    layout.addRead(reads[readIndex], minBaseQuality)
                }
                // Next iteration, we may be able to overlap more reads now.
            }
        }
    }

    companion object
    {
        private val sLogger = LogManager.getLogger(VDJSequenceBuilder::class.java)

        private fun mergeIdentical(vjGeneToVdjMap: Map<Pair<VJGeneType?, VJGeneType?>, List<VDJSequence>>, minBaseQuality: Byte,
                                   threadCount: Int): List<VDJSequence>
        {
            sLogger.debug("Merging identical VDJ sequences")

            // Naive multithreading because otherwise this step can take over an hour on samples with many sequences.
            // It's likely that a few genes have most of the sequences, which will limit the parallel speedup, but it's better than nothing.

            // Note we are careful to collect the results in a fixed order to produce deterministic results.
            val threadResults = vjGeneToVdjMap.entries.map { ArrayList<VDJSequence>() }.toMutableList()
            val runnables = vjGeneToVdjMap.entries.sortedWith(compareBy({ it.key.first }, { it.key.second }))
                .map { it.value }
                .withIndex()
                .map { (index, vdjs) ->
                    Runnable {
                        val merged = mergeIdentical(vdjs, minBaseQuality)
                        synchronized(threadResults) {
                            threadResults[index] = merged
                        }
                    }}
            TaskExecutor.executeRunnables(runnables, threadCount)
            return threadResults.flatten()
        }

        // for VDJ sequences that are identical, we put them together
        private fun mergeIdentical(vdjList: List<VDJSequence>, minBaseQuality: Byte) : ArrayList<VDJSequence>
        {
            /*
            // first parse we do it in a fast way to find the ones that are exactly the same
            val seqHashMap: ArrayListMultimap<String, VDJSequence> = ArrayListMultimap.create()

            vdjList.forEach({ vdj -> seqHashMap.put(vdj.sequence, vdj) })

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
                val diffAccumulator = { s1: Map.Entry<Byte, Int>, s2: Map.Entry<Byte, Int> ->
                    if (s1.value >= 1 && s2.value >= 1)
                    {
                        ++numBaseDiff
                        maxReadsAtDiff = max(maxReadsAtDiff, min(s1.value, s2.value))
                    }

                    // we say they are the same if there is only 1 base difference
                    // and only 1 read supporting that one base diff
                    (numBaseDiff <= 1 && maxReadsAtDiff < 2)
                }

                if (VdjBuilderUtils.vdjSequenceCompare(vdj1, vdj2, diffAccumulator))
                {
                    sLogger.trace("removing {} as it is similar to {}",
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

        private fun vdjSequenceIdentical(vdj1: VDJSequence, vdj2: VDJSequence) : Boolean
        {
            var numBaseDiff = 0
            val diffAccumulator = { baseSupport1: Map.Entry<Byte, Int>, baseSupport2: Map.Entry<Byte, Int> ->

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