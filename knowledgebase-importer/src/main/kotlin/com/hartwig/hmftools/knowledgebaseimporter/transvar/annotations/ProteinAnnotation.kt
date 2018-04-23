package com.hartwig.hmftools.knowledgebaseimporter.transvar.annotations

data class ProteinAnnotation(override val transcript: String, override val alteration: String) : Annotation
