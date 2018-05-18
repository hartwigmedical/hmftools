package com.hartwig.hmftools.knowledgebaseimporter.output

import com.hartwig.hmftools.common.variant.SomaticVariant

data class SomaticVariantEvent(val chromosome: String, val position: String, val ref: String, val alt: String) :
        ActionableEvent {
    companion object {
        val header = listOf("chromosome", "position", "ref", "alt")

        operator fun invoke(variant: SomaticVariant): SomaticVariantEvent {
            return SomaticVariantEvent(variant.chromosome(), variant.position().toString(), variant.ref(), variant.alt())
        }
    }

    val record: List<String> = listOf(chromosome, position, ref, alt)

    override fun toString(): String {
        return "$chromosome:$position $ref->$alt"
    }
}
