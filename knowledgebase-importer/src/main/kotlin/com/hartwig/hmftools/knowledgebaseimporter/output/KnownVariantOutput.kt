package com.hartwig.hmftools.knowledgebaseimporter.output

import com.hartwig.hmftools.common.variant.SomaticVariant

data class KnownVariantOutput(private val gene: String, private val transcript: String, private val additionalInfo: String,
                              val variant: SomaticVariant) {
    enum class Header {
        GENE, TRANSCRIPT, INFO, CHROMOSOME, POSITION, REF, ALT
    }

    val record: List<String> = listOf(gene,
                                      transcript,
                                      additionalInfo,
                                      variant.chromosome(),
                                      variant.position().toString(),
                                      variant.ref(),
                                      variant.alt())
}
