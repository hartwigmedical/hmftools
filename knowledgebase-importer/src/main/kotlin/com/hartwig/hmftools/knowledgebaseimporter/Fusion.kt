package com.hartwig.hmftools.knowledgebaseimporter

sealed class Fusion
data class FusionPair(val fiveGene: String, val threeGene: String) : Fusion()
data class PromiscuousGene(val gene: String) : Fusion()
