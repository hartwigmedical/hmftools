package com.hartwig.hmftools.knowledgebaseimporter.cgi.input

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.KnowledgebaseEvent

interface CgiInput : KnowledgebaseEvent {
    val gdna: String
}
