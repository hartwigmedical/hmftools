package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordMetadata

data class OncoKbMetadata(override val gene: String, override val transcript: String) : RecordMetadata {
    override val source = "oncoKb"
}
