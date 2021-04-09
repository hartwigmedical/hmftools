package com.hartwig.hmftools.lilac.variant

import com.hartwig.hmftools.common.variant.CodingEffect
import com.hartwig.hmftools.lilac.hla.HlaAllele

data class SomaticCodingCount(val allele: HlaAllele, val missense: Double, val nonsense: Double, val splice: Double, val synonmous: Double) {
    val total = missense + nonsense + splice + synonmous

    companion object {
        fun create(winners: List<HlaAllele>): List<SomaticCodingCount> {
            val result = winners.map { SomaticCodingCount(it, 0.0, 0.0, 0.0, 0.0) }
            return result
        }

        fun List<SomaticCodingCount>.addVariant(effect: CodingEffect, alleles: List<HlaAllele>): List<SomaticCodingCount> {
            return this.map { it.addVariant(effect, alleles) }
        }
    }

    fun addVariant(effect: CodingEffect, alleles: List<HlaAllele>): SomaticCodingCount {
        if (allele !in alleles) {
            return this
        }

        val contribution = 1.0 / alleles.size
        return when (effect) {
            CodingEffect.MISSENSE -> copy(missense = missense + contribution)
            CodingEffect.NONSENSE_OR_FRAMESHIFT -> copy(nonsense = nonsense + contribution)
            CodingEffect.SPLICE -> copy(splice = splice + contribution)
            CodingEffect.SYNONYMOUS -> copy(synonmous = synonmous + contribution)
            else -> this
        }
    }

}