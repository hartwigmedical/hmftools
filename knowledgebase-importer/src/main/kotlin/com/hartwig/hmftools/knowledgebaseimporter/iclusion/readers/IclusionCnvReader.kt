package com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers

import com.hartwig.hmftools.knowledgebaseimporter.iclusion.IclusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.EventType
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader
import com.hartwig.hmftools.knowledgebaseimporter.output.CnvEvent

object IclusionCnvReader : SomaticEventReader<IclusionEvent, CnvEvent> {
    private fun match(event: IclusionEvent): Boolean {
        return event.types.size == 1 && event.types.contains(EventType.CNV)
    }

    override fun read(event: IclusionEvent): List<CnvEvent> {
        if (match(event)) return listOfNotNull(readEvent(event))
        return emptyList()
    }

    private fun readEvent(event: IclusionEvent): CnvEvent? {
        return when (event.variant) {
            "AMPLIFICATION" -> CnvEvent.amplification(event.gene)
            "COPY-GAIN"     -> CnvEvent.amplification(event.gene)
            "DELETION"      -> CnvEvent.deletion(event.gene)
            "LOSS"          -> CnvEvent.deletion(event.gene)
            else            -> null
        }
    }
}
