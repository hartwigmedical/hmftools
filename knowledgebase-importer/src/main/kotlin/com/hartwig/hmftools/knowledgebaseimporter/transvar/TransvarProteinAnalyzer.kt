package com.hartwig.hmftools.knowledgebaseimporter.transvar

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ProteinAnnotation

data class TransvarProteinAnalyzer(override val transvarLocation: String) : TransvarAnalyzer<ProteinAnnotation> {
    override val analysisMode: String = "panno"
}
