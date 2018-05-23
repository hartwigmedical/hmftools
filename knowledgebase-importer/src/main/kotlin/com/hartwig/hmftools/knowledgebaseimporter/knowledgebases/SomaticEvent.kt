package com.hartwig.hmftools.knowledgebaseimporter.knowledgebases

interface SomaticEvent

data class GDnaVariant(val gDnaImpact: String) : SomaticEvent
data class KnowledgebaseVariant(val chromosome: String, val position: Long, val ref: String?, val alt: String?) : SomaticEvent
