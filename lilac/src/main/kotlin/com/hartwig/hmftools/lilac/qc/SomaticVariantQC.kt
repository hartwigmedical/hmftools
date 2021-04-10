package com.hartwig.hmftools.lilac.qc

import com.hartwig.hmftools.lilac.variant.SomaticCodingCount

data class SomaticVariantQC(val variantCount: Int, val variantAlleleCount: Double) {

    fun header(): List<String> {
        return listOf("variantCount", "variantAlleleCount")
    }

    fun body(): List<String> {
        return listOf(variantCount.toString(), variantAlleleCount.toString())
    }

    companion object {
        fun create(variantCount: Int, codingCount: List<SomaticCodingCount>): SomaticVariantQC {
            val totalCount = codingCount.map { it.total }.sum()
            return SomaticVariantQC(variantCount, totalCount)
        }
    }

}