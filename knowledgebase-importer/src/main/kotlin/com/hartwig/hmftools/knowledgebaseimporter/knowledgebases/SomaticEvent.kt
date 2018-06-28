package com.hartwig.hmftools.knowledgebaseimporter.knowledgebases

interface SomaticEvent

data class GDnaVariant(val gDnaImpact: String) : SomaticEvent

data class KnowledgebaseVariant(val gene: String, val chromosome: String, val position: Long, val ref: String?, val alt: String?) :
        SomaticEvent

interface GenericMutation : SomaticEvent {
    val transcript: String?
    val gene: String
}

sealed class RangeMutation : GenericMutation

data class CodonMutations(override val gene: String, override val transcript: String?, val codonNumber: Int) : GenericMutation,
        RangeMutation()

data class CodonRangeMutations(override val gene: String, override val transcript: String?, val startCodon: Int,
                               val endCodon: Int) : GenericMutation, RangeMutation()

data class GenericRangeMutations(override val gene: String, override val transcript: String?, val startPosition: Int,
                                 val endPosition: Int) : GenericMutation, RangeMutation()

data class ExonMutations(override val gene: String, override val transcript: String?, val exonNumber: Int) : GenericMutation,
        RangeMutation()

sealed class HgvsAnnotation : SomaticEvent {
    abstract val transcript: String
    abstract val alteration: String
}

data class CDnaAnnotation(override val transcript: String, override val alteration: String) : HgvsAnnotation()

data class ProteinAnnotation(override val transcript: String, override val alteration: String) : HgvsAnnotation()
