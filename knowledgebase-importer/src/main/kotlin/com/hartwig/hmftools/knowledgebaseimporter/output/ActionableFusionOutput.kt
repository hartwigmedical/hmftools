package com.hartwig.hmftools.knowledgebaseimporter.output

data class ActionableFusionOutput(val fusion: Fusion, val actionability: Actionability) {
    companion object {
        val fusionPairHeader = listOf("fiveGene", "threeGene") + Actionability.header
        val promiscuousGeneHeader = listOf("gene") + Actionability.header
    }

    val record: List<String> = fusion.record + actionability.record
}
