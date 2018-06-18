package com.hartwig.hmftools.knowledgebaseimporter.knowledgebases

interface SomaticEvent

data class GDnaVariant(val gDnaImpact: String) : SomaticEvent

data class KnowledgebaseVariant(val gene: String, val chromosome: String, val position: Long, val ref: String?, val alt: String?) :
        SomaticEvent

sealed class GenericMutation : SomaticEvent {
    abstract val transcript: String?
    abstract val gene: String
}

data class CodonMutations(override val gene: String, override val transcript: String?, val codonNumber: Int) : GenericMutation()

data class CodonRangeMutations(override val gene: String, override val transcript: String?, val startCodon: Int,
                               val endCodon: Int) : GenericMutation()

data class ExonMutations(override val gene: String, override val transcript: String?, val exonNumber: Int) : GenericMutation()

data class GeneMutations(override val gene: String, override val transcript: String?) : GenericMutation()

sealed class HgvsAnnotation : SomaticEvent {
    abstract val transcript: String
    abstract val alteration: String
}

data class CDnaAnnotation(override val transcript: String, override val alteration: String) : HgvsAnnotation()

data class ProteinAnnotation(override val transcript: String, override val alteration: String) : HgvsAnnotation()
