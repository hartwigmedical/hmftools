package com.hartwig.hmftools.knowledgebaseimporter.output

sealed class Fusion {
    abstract val record: List<String>
}

data class FusionPair(val fiveGene: String, val threeGene: String) : Fusion() {
    companion object {
        val header = listOf("fiveGene", "threeGene")
    }

    override val record: List<String> = listOf(fiveGene, threeGene)
}

data class PromiscuousGene(val gene: String) : Fusion() {
    companion object {
        val header = listOf("gene")
    }

    override val record: List<String> = listOf(gene)
}
