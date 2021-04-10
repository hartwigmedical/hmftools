package com.hartwig.hmftools.lilac.qc

import com.hartwig.hmftools.lilac.variant.SomaticCodingCount
import org.apache.logging.log4j.LogManager

data class SomaticVariantQC(val variantCount: Int, val variantAlleleCount: Double) {

    fun header(): List<String> {
        return listOf("variantCount", "variantAlleleCount")
    }

    fun body(): List<String> {
        return listOf(variantCount.toString(), variantAlleleCount.toString())
    }

    companion object {
        val logger = LogManager.getLogger(this::class.java)

        fun create(variantCount: Int, codingCount: List<SomaticCodingCount>): SomaticVariantQC {
            val totalCount = codingCount.map { it.total }.sum()
            if (Math.abs(totalCount - variantCount) > 0.01) {
                logger.warn("    UNASSIGNED_VARIANT - $variantCount variants found but $totalCount assigned")
            }

            return SomaticVariantQC(variantCount, totalCount)
        }
    }

}