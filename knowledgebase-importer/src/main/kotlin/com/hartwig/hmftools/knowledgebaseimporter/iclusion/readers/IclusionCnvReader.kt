package com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers

import com.hartwig.hmftools.knowledgebaseimporter.iclusion.IclusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader
import com.hartwig.hmftools.knowledgebaseimporter.output.CnvEvent

object IclusionCnvReader : SomaticEventReader<IclusionEvent, CnvEvent> {

    override fun read(event: IclusionEvent): List<CnvEvent> {
        return listOfNotNull(readEvent(event))
    }

    private fun readEvent(event: IclusionEvent): CnvEvent? {
        return when (event.variant) {
            "AMPLIFICATION"     -> CnvEvent.amplification(event.gene)
            "COPY-GAIN"         -> CnvEvent.amplification(event.gene)
            "OVEREXPRESSION"    -> CnvEvent.amplification(event.gene)
            "DELETION"          -> CnvEvent.deletion(event.gene)
            "LOSS"              -> CnvEvent.deletion(event.gene)
            "LOSS-OF-FUNCTION"  -> CnvEvent.deletion(event.gene)
            else                -> null
        }
    }
}
