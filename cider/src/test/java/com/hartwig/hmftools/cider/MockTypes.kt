package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.genes.VJ
import com.hartwig.hmftools.cider.genes.VJGeneType
import com.hartwig.hmftools.cider.layout.ReadLayout
import java.util.*

// create a mock anchor blosum searcher that returns the match that we give it
class MockAnchorBlosumSearcher : IAnchorBlosumSearcher
{
    var anchorBlosumMatch: AnchorBlosumMatch? = null

    override fun searchForAnchor(readString: String, mode: IAnchorBlosumSearcher.Mode) : AnchorBlosumMatch?
    {
        return anchorBlosumMatch
    }

    override fun searchForAnchor(sequence: String, targetAnchorGeneTypes: Collection<VJGeneType>,
                                 mode: IAnchorBlosumSearcher.Mode,
                                 startOffset: Int, endOffset: Int) : AnchorBlosumMatch?
    {
        return anchorBlosumMatch
    }
}

// create a mock layout adaptor that just return anchor range that we give it
class MockVJReadLayoutAdaptor : IVJReadLayoutAdaptor()
{
    val anchorRangeMap = IdentityHashMap<ReadLayout, IntRange>()
    val readCandidateMap = IdentityHashMap<ReadLayout.Read, VJReadCandidate>()
    val readSliceMap = IdentityHashMap<ReadLayout.Read, ReadSlice>()
    // val templateSequenceMap = IdentityHashMap<ReadLayout.Read, String>()

    override fun toReadCandidate(read: ReadLayout.Read) : VJReadCandidate
    {
        return readCandidateMap[read]!!
    }

    override fun toLayoutReadSlice(read: ReadLayout.Read) : ReadSlice
    {
        return readSliceMap[read]!!
    }

    override fun getAnchorMatchMethod(layout: ReadLayout): VJReadCandidate.MatchMethod
    {
        return VJReadCandidate.MatchMethod.ALIGN
    }

    override fun getTemplateAnchorSequence(layout: ReadLayout) : String
    {
        return "TGACCC"
    }

    override fun getExtrapolatedAnchorRange(vj: VJ, layout: ReadLayout) : IntRange
    {
        return anchorRangeMap[layout]!!
    }
}
