package com.hartwig.hmftools.cider.primer

import com.hartwig.hmftools.cider.VJReadLayoutAdaptor
import com.hartwig.hmftools.cider.VDJSequence
import com.hartwig.hmftools.cider.layout.ReadLayout
import org.apache.logging.log4j.LogManager

// match the VDJ sequence with the primers
class VdjPrimerMatcher(private val vjLayoutAdaptor: VJReadLayoutAdaptor, private val maxMismatch: Int)
{
    fun matchVdjPrimer(vdjList: List<VDJSequence>, primerList: List<Primer>) : List<VdjPrimerMatch>
    {
        val matchList = ArrayList<VdjPrimerMatch>()

        for (vdj in vdjList)
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

                    matchList.add(VdjPrimerMatch(vdj, primer, i, sequence))
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
