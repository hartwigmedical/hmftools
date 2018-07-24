package com.hartwig.hmftools.knowledgebaseimporter.cgi

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordMetadata

data class CgiMetadata(override val gene: String, override val transcript: String) : RecordMetadata {
    override val source = "cgi"
}
