package com.hartwig.hmftools.knowledgebaseimporter.output

sealed class FusionEvent : ActionableEvent {
    abstract val record: List<String>
}

data class FusionPair(val fiveGene: String, val threeGene: String) : FusionEvent() {
    companion object {
        val header = listOf("fiveGene", "threeGene")
    }

    override val record: List<String> = listOf(fiveGene, threeGene)

    override fun toString(): String {
        return "$fiveGene - $threeGene fusion"
    }
}

data class PromiscuousGene(val gene: String) : FusionEvent() {
    companion object {
        val header = listOf("gene")
    }

    override val record: List<String> = listOf(gene)

    override fun toString(): String {
        return "$gene fusions"
    }
}
