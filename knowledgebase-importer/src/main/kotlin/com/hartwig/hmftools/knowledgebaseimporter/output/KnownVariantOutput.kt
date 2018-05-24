package com.hartwig.hmftools.knowledgebaseimporter.output

data class KnownVariantOutput(private val transcript: String, private val additionalInfo: String,
                              private val variant: SomaticVariantEvent) {
    companion object {
        val header = listOf("gene", "transcript", "info") + SomaticVariantEvent.header
    }

    val record: List<String> = listOf(variant.gene, transcript, additionalInfo) + variant.record
}
