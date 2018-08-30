package com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ProteinAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.KnowledgebaseEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader
import com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers.TransvarProteinMatcher

object IclusionProteinAnnotationReader : SomaticEventReader<KnowledgebaseEvent, ProteinAnnotation> {
    private fun match(event: KnowledgebaseEvent) = TransvarProteinMatcher.matches(event.variant)

    override fun read(event: KnowledgebaseEvent): List<ProteinAnnotation> {
        val transcript = event.transcript
        if (match(event) && transcript != null)
            return listOfNotNull(ProteinAnnotation(transcript, event.variant, TransvarProteinMatcher.type(event.variant)))
        return emptyList()
    }
}
