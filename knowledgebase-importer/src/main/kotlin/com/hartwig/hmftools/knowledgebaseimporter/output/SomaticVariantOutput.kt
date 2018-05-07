package com.hartwig.hmftools.knowledgebaseimporter.output

import com.hartwig.hmftools.common.variant.SomaticVariant

data class SomaticVariantOutput(private val chromosome: String, private val position: String, private val ref: String,
                                private val alt: String) {
    companion object {
        val header = listOf("chromosome", "position", "ref", "alt")

        operator fun invoke(variant: SomaticVariant): SomaticVariantOutput {
            return SomaticVariantOutput(variant.chromosome(), variant.position().toString(), variant.ref(), variant.alt())
        }
    }

    val record: List<String> = listOf(chromosome, position, ref, alt)
}
