package com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers

import com.hartwig.hmftools.knowledgebaseimporter.FusionReader
import com.hartwig.hmftools.knowledgebaseimporter.iclusion.IclusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.EventType
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionEvent

object IclusionFusionReader : SomaticEventReader<IclusionEvent, FusionEvent> {
    const val FUSION_SEPARATOR = "-"
    private val fusionReader = FusionReader(separators = setOf(FUSION_SEPARATOR))

    private fun match(event: IclusionEvent): Boolean {
        return (event.types.size == 1 && event.types.contains(EventType.FUS))
    }

    override fun read(event: IclusionEvent): List<FusionEvent> {
        if (match(event)) return listOfNotNull(fusionReader.read(event.gene, event.variant))
        return emptyList()
    }
}
