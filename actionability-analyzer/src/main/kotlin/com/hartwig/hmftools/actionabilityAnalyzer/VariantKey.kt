package com.hartwig.hmftools.actionabilityAnalyzer

import com.hartwig.hmftools.knowledgebaseimporter.output.SomaticVariantEvent
import com.hartwig.hmftools.patientdb.data.PotentialActionableVariant

data class VariantKey(val chromosome: String, val position: String, val ref: String, val alt: String) {
    companion object {
        operator fun invoke(variant: PotentialActionableVariant): VariantKey {
            return VariantKey(variant.chromosome(), variant.position().toString(), variant.ref(), variant.alt())
        }

        operator fun invoke(variant: SomaticVariantEvent): VariantKey {
            return VariantKey(variant.chromosome, variant.position, variant.ref, variant.alt)
        }
    }
}
