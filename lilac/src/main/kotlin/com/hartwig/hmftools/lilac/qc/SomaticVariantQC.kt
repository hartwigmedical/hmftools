package com.hartwig.hmftools.lilac.qc

import com.hartwig.hmftools.lilac.variant.SomaticCodingCount
import org.apache.logging.log4j.LogManager
import kotlin.math.abs

data class SomaticVariantQC(val variantCount: Int, val variantAlleleCount: Double) {

    fun header(): List<String> {
        return listOf("variantCount", "variantAlleleCount")
    }

    fun body(): List<String> {
        return listOf(variantCount.toString(), variantAlleleCount.toString())
    }

    fun unmatchedVariants(): Boolean {
        return abs(variantCount - variantCount) > 0.01
    }

    companion object {
        val logger = LogManager.getLogger(this::class.java)

        fun create(variantCount: Int, codingCount: List<SomaticCodingCount>): SomaticVariantQC {
            val totalCount = codingCount.map { it.total }.sum()
            val result =  SomaticVariantQC(variantCount, totalCount)
            if (result.unmatchedVariants()) {
                logger.warn("    UNASSIGNED_VARIANT - $variantCount variants found but $totalCount assigned")
            }

            return result;
        }
    }

}