package com.hartwig.hmftools.knowledgebaseimporter.knowledgebases

interface SomaticEvent

data class GDnaVariant(val gDnaImpact: String) : SomaticEvent

data class KnowledgebaseVariant(val gene: String, val chromosome: String, val position: Long, val ref: String?, val alt: String?) :
        SomaticEvent

sealed class GenericMutation : SomaticEvent

data class CodonMutations(val gene: String, val transcript: String?, val codonNumber: Int) : GenericMutation()

data class CodonRangeMutations(val gene: String, val transcript: String?, val startCodon: Int, val endCodon: Int) : GenericMutation()

data class ExonMutations(val gene: String, val transcript: String?, val exonNumber: Int) : GenericMutation()

data class GeneMutations(val gene: String, val transcript: String?) : GenericMutation()

sealed class HgvsAnnotation : SomaticEvent {
    abstract val transcript: String
    abstract val alteration: String
}

data class CDnaAnnotation(override val transcript: String, override val alteration: String) : HgvsAnnotation()

data class ProteinAnnotation(override val transcript: String, override val alteration: String) : HgvsAnnotation()
