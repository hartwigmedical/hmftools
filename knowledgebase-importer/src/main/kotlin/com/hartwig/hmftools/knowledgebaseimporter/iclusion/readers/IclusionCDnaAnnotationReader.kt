package com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CDnaAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.KnowledgebaseEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader
import com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers.TransvarCDnaMatcher

object IclusionCDnaAnnotationReader : SomaticEventReader<KnowledgebaseEvent, CDnaAnnotation> {
    private fun match(event: KnowledgebaseEvent) = TransvarCDnaMatcher.matches(event.variant)

    override fun read(event: KnowledgebaseEvent): List<CDnaAnnotation> {
        val transcript = event.transcript
        if (match(event) && transcript != null)
            return listOfNotNull(CDnaAnnotation(transcript, event.variant, TransvarCDnaMatcher.type(event.variant)))
        return emptyList()
    }
}
