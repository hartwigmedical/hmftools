package com.hartwig.hmftools.knowledgebaseimporter.transvar

import com.hartwig.hmftools.knowledgebaseimporter.transvar.annotations.CDnaAnnotation

data class TransvarCdnaAnalyzer(override val transvarLocation: String) : TransvarAnalyzer<CDnaAnnotation> {
    override val analysisMode: String = "canno"
}
