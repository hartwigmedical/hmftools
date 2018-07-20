package com.hartwig.hmftools.actionabilityAnalyzer

import com.hartwig.hmftools.knowledgebaseimporter.output.SomaticVariantEvent

data class VariantKey(val chromosome: String, val position: String, val ref: String, val alt: String) {
    companion object {
        operator fun invoke(variant: SomaticVariantEvent): VariantKey {
            return VariantKey(variant.chromosome, variant.position, variant.ref, variant.alt)
        }

        operator fun invoke(variant: CohortMutation): VariantKey {
            return VariantKey(variant.chromosome, variant.position, variant.ref, variant.alt)
        }
    }
}
