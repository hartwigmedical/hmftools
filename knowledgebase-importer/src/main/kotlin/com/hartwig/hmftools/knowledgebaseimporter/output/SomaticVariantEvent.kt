package com.hartwig.hmftools.knowledgebaseimporter.output

import com.hartwig.hmftools.common.variant.SomaticVariant

data class SomaticVariantEvent(val gene: String, val chromosome: String, val position: String, val ref: String, val alt: String) :
        ActionableEvent {
    companion object {
        val header = listOf("chromosome", "position", "ref", "alt")

        operator fun invoke(gene: String, variant: SomaticVariant): SomaticVariantEvent {
            return SomaticVariantEvent(gene, variant.chromosome(), variant.position().toString(), variant.ref(), variant.alt())
        }
    }

    val record: List<String> = listOf(chromosome, position, ref, alt)

    override fun toString(): String {
        return "$gene $chromosome:$position $ref->$alt"
    }
}
