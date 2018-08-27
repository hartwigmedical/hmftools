package com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers

import com.hartwig.hmftools.knowledgebaseimporter.iclusion.IclusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CDnaAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.EventType
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader
import com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers.TransvarCDnaMatcher

object IclusionCDnaAnnotationReader : SomaticEventReader<IclusionEvent, CDnaAnnotation> {
    private fun match(event: IclusionEvent): Boolean {
        return event.types.size == 1 && event.types.contains(EventType.MUT) && TransvarCDnaMatcher.matches(event.variant)
    }

    override fun read(event: IclusionEvent): List<CDnaAnnotation> {
        if (match(event)) return listOfNotNull(CDnaAnnotation(event.transcript, event.variant))
        return emptyList()
    }
}
