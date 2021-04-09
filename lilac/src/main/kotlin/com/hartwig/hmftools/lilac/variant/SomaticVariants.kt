package com.hartwig.hmftools.lilac.variant

import com.hartwig.hmftools.common.variant.CodingEffect
import com.hartwig.hmftools.common.variant.VariantContextDecorator
import com.hartwig.hmftools.lilac.LilacApplication
import com.hartwig.hmftools.lilac.LilacConfig
import htsjdk.variant.vcf.VCFFileReader
import org.apache.logging.log4j.LogManager
import java.io.File

class SomaticVariants(private val config: LilacConfig) {

    private val logger = LogManager.getLogger(this::class.java)
    private val UNKNOWN_CODING_EFFECT = setOf(CodingEffect.NONE, CodingEffect.UNDEFINED)

    fun readSomaticVariants(): List<VariantContextDecorator> {
        if (config.tumorBam.isEmpty() || config.somaticVcf.isEmpty()) {
            return listOf()
        }

        val contig = LilacApplication.HLA_TRANSCRIPTS.map { it.chromosome() }.first()
        val minPosition = LilacApplication.HLA_TRANSCRIPTS.map { it.start() }.min()!!
        val maxPosition = LilacApplication.HLA_TRANSCRIPTS.map { it.end() }.max()!!

        logger.info("Reading somatic vcf: ${config.somaticVcf}")
        val result = mutableListOf<VariantContextDecorator>()
        val fileReader = VCFFileReader(File(config.somaticVcf), false)
        val iterator = if (fileReader.isQueryable) {
            fileReader.query(contig, minPosition.toInt(), maxPosition.toInt())
        } else (fileReader.iterator())

        for (variantContext in iterator) {
            val enriched = VariantContextDecorator(variantContext)
            if (enriched.gene() in LilacApplication.HLA_GENES && enriched.canonicalCodingEffect() !in UNKNOWN_CODING_EFFECT && enriched.isPass) {
                result.add(enriched)
            }
        }

        fileReader.close()

        logger.info("    found ${result.size} HLA variants")
        return result
    }

}