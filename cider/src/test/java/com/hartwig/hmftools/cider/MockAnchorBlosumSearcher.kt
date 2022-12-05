package com.hartwig.hmftools.cider

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
