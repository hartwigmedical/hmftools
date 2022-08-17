package com.hartwig.hmftools.cider

// create a mock anchor blosum searcher that returns the match that we give it
class MockAnchorBlosumSearcher : IAnchorBlosumSearcher
{
    var anchorBlosumMatch: AnchorBlosumMatch? = null

    override fun searchForAnchor(dnaSeq: String, targetAnchorGeneType: VJGeneType, startOffset: Int,
                                 endOffset: Int) : AnchorBlosumMatch?
    {
        return anchorBlosumMatch
    }
}
