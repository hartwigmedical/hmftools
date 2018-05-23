package com.hartwig.hmftools.knowledgebaseimporter.knowledgebases

import com.hartwig.hmftools.knowledgebaseimporter.output.Actionability

interface ActionableRecord : KnowledgebaseRecord {
    val actionability: List<Actionability>
}
