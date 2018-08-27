package com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events

interface KnowledgebaseEvent {
    val gene: String
    val variant: String
    val transcript: String?
}
