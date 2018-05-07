package com.hartwig.hmftools.knowledgebaseimporter.output

data class KnownVariantOutput(private val gene: String, private val transcript: String, private val additionalInfo: String,
                              private val variant: SomaticVariantOutput) {
    companion object {
        val header = listOf("gene", "transcript", "info") + SomaticVariantOutput.header
    }

    val record: List<String> = listOf(gene, transcript, additionalInfo) + variant.record
}
