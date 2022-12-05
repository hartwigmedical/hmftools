package com.hartwig.hmftools.cider.primer

import com.hartwig.hmftools.cider.VDJSequence
import com.hartwig.hmftools.cider.layout.ReadLayout
import org.apache.logging.log4j.LogManager

// match the VDJ sequence with the primers
class VdjPrimerMatcher(private val maxMismatch: Int)
{
    fun matchVdjPrimer(vdjList: List<VDJSequence>, primerList: List<Primer>) : List<VdjPrimerMatch>
    {
        val matchList = ArrayList<VdjPrimerMatch>()

        for (vdj in vdjList)
        {
            // in order to match primer list, we need to build up the layout again
            val fullLayout: ReadLayout = vdj.layout

            // now go through each list and see if anything matches
            val sequence = fullLayout.consensusSequence()

            for (primer in primerList)
            {
                val primerSeq = primer.sequence

                for (i in 0 .. sequence.length - primerSeq.length)
                {
                    // try to find it using a simple loop
                    var numMismatch = 0

                    for (j in primerSeq.indices)
                    {
                        if (sequence[i + j] != primerSeq[j])
                        {
                            ++numMismatch

                            if (numMismatch > maxMismatch)
                            {
                                break
                            }
                        }
                    }

                    if (numMismatch <= maxMismatch)
                    {
                        sLogger.info("vdj seq matches primer: i: {}, {}, full vdj seq: {}, mismatch: {}",
                            i, primer, vdj.sequenceFormatted, numMismatch)

                        matchList.add(VdjPrimerMatch(vdj, primer, i, numMismatch, sequence))
                        break
                    }
                }
            }
        }

        return matchList
    }

    /*
    fun matchPrimerList(vdj: VDJSequence, primerList: List<Primer>)
    {
        // in order to match primer list, we need to build up the layout again
        val fullLayout: ReadLayout = vjLayoutAdaptor.buildFullLayout(vdj.layout)

        // now go through each list and see if anything matches
        val sequence = fullLayout.consensusSequence()

        for (primer in primerList)
        {
            val i = sequence.indexOf(primer.sequence)
            if (i != -1)
            {
                sLogger.info("vdj seq matches primer: {}, full vdj seq: {}", primer, vdj.sequenceFormatted)
            }
        }
    }*/

    companion object
    {
        private val sLogger = LogManager.getLogger(VdjPrimerMatcher::class.java)
    }
}
