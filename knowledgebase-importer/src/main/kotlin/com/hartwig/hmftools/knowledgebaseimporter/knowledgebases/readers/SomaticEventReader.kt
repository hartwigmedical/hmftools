package com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.KnowledgebaseEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent

interface SomaticEventReader<in R : KnowledgebaseEvent, out E : SomaticEvent> {
    fun read(event: R): List<E>
}
