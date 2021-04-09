package com.hartwig.hmftools.lilac.qc

import com.hartwig.hmftools.lilac.variant.SomaticCodingCount

data class SomaticVariantQC(val variantCount: Int, val variantAlleleCount: Double) {

    companion object {
        fun create(variantCount: Int, codingCount: List<SomaticCodingCount>) {
            val totalCoverage = codingCount.map { it.total }.sum()

        }

    }

}