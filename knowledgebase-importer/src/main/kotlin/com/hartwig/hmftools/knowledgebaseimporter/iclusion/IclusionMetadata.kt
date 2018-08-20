package com.hartwig.hmftools.knowledgebaseimporter.iclusion

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordMetadata

data class IclusionMetadata(override val gene: String, override val transcript: String) : RecordMetadata {
    override val source = "iclusion"
}
