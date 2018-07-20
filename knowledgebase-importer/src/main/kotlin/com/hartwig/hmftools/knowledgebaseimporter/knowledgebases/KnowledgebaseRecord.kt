package com.hartwig.hmftools.knowledgebaseimporter.knowledgebases

interface KnowledgebaseRecord : RecordMetadata {
    val events: List<SomaticEvent>
}
