package com.hartwig.hmftools.cider

/*
// try to overlay V and J together and create a single VDJ
class VJLayoutAligner(minMatchedBases: Int, minMatchRatio: Double)
{
    val minMatchedBases: Int = minMatchedBases
    val minMatchRatio: Double = minMatchRatio

    fun findAlignedVJs(vLayouts: List<ReadLayout>, jLayouts: List<ReadLayout>): List<VDJCandidate>
    {
        val alignedVJs = ArrayList<VDJCandidate>()

        // we try to align each pair together
        for (vLayout in vLayouts)
        {
            for (jLayout in jLayouts)
            {
                val vdjCandidate: VDJCandidate? = tryAlignVJ(vLayout, jLayout)

                if (vdjCandidate != null)
                    alignedVJs.add(vdjCandidate)
            }
        }

        return alignedVJs
    }

    // we should probably say allow some more leeway if a position has very little support
    fun tryAlignVJ(vLayout: ReadLayout, jLayout: ReadLayout): VDJCandidate?
    {
        val vSupport = vLayout.highQualSequenceSupport.support
        val jSupport = jLayout.highQualSequenceSupport.support

        // try multiple location to align them
        for (i in minMatchedBases until Math.min(vSupport.size, jSupport.size))
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
            }
        }

        /*val fullSeq = vSeq.substring(0, alignment.seq1AlignStart) + overlapSeqBuilder.toString() + jSeq.substring(alignment.seq2AlignEnd)

        sLogger.debug("v and j overlap")
        sLogger.debug("v: {}", vSeq)
        sLogger.debug("j: {}", jSeq)
        sLogger.debug("overlap: v({}), j({})", vOverlapSeq, jOverlapSeq)
        sLogger.debug("full seq: {}", fullSeq)

        return VDJCandidate(vLayout, jLayout, alignment.seq1AlignStart, alignment.seq2AlignEnd, overlapSeqBuilder.toString())*/
        return null
    }

    companion object
    {
        private val sLogger = LogManager.getLogger(VJLayoutAligner::class.java)
    }
}
 */