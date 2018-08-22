package com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.KnowledgebaseEvent

interface SomaticEventReader<in R : KnowledgebaseEvent, out E : SomaticEvent> {
    fun read(event: R): List<E>
}
