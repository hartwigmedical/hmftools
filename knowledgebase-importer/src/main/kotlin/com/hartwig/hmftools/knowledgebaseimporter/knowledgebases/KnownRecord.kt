package com.hartwig.hmftools.knowledgebaseimporter.knowledgebases

interface KnownRecord : KnowledgebaseRecord {
    val reference: String
    val annotation: String
}
