package com.hartwig.hmftools.knowledgebaseimporter.knowledgebases

interface KnowledgebaseEvent {
    val source: String
    val gene: String
    val transcript: String
    val variant: String
    val types: List<EventType>
}
