package com.hartwig.hmftools.knowledgebaseimporter.transvar.annotations

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent

interface Annotation : SomaticEvent {
    val transcript: String
    val alteration: String
}
