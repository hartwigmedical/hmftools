package com.hartwig.hmftools.knowledgebaseimporter.knowledgebases

interface SomaticEvent

sealed class HgvsAnnotation : SomaticEvent {
    abstract val transcript: String
    abstract val alteration: String
}

data class GDnaVariant(val gDnaImpact: String) : SomaticEvent

data class KnowledgebaseVariant(val gene: String, val chromosome: String, val position: Long, val ref: String?, val alt: String?) :
        SomaticEvent

data class CDnaAnnotation(override val transcript: String, override val alteration: String) : HgvsAnnotation()

data class ProteinAnnotation(override val transcript: String, override val alteration: String) : HgvsAnnotation()
