package com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers

import com.hartwig.hmftools.knowledgebaseimporter.iclusion.IclusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.EventType
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
            "AMPLIFICATION" -> CnvEvent(event.gene, "Amplification")
            "COPY-GAIN"     -> CnvEvent(event.gene, "Amplification")
            "DELETION"      -> CnvEvent(event.gene, "Deletion")
            "LOSS"          -> CnvEvent(event.gene, "Deletion")
            else            -> null
        }
    }
}
