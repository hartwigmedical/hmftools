package com.hartwig.hmftools.knowledgebaseimporter.civic

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordMetadata

data class CivicMetadata(override val gene: String, override val transcript: String) : RecordMetadata {
    override val source = "civic"
}
