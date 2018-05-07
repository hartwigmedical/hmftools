package com.hartwig.hmftools.knowledgebaseimporter.output

data class ActionableVariantOutput(private val gene: String, private val variant: SomaticVariantOutput,
                                   private val actionability: Actionability) {
    companion object {
        val header = listOf("gene") + SomaticVariantOutput.header + Actionability.header
    }

    val record: List<String> = listOf(gene) + variant.record + actionability.record
}
