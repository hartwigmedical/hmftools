package com.hartwig.hmftools.lilac.variant

import com.hartwig.hmftools.common.variant.CodingEffect
import com.hartwig.hmftools.common.variant.VariantContextDecorator
import com.hartwig.hmftools.lilac.LilacApplication
import com.hartwig.hmftools.lilac.LilacConfig
import com.hartwig.hmftools.lilac.LociPosition
import com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage
import com.hartwig.hmftools.lilac.read.FragmentAlleles
import com.hartwig.hmftools.lilac.read.SAMRecordReader
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.vcf.VCFFileReader
import org.apache.logging.log4j.LogManager
import java.io.File

class Somatics {

    private val logger = LogManager.getLogger(this::class.java)

    companion object {
        val UNKNOWN_CODING_EFFECT = setOf(CodingEffect.NONE, CodingEffect.UNDEFINED)
    }


    fun process(config: LilacConfig, reader: SAMRecordReader, winners: List<HlaSequenceLoci>, hetLoci: Collection<Int>, lociPosition: LociPosition) {
        if (config.tumorBam.isEmpty() || config.somaticVcf.isEmpty()) {
            return
        }

        val variants = readSomatics(config)
        if (variants.isEmpty()) {
            return
        }

        val variantLoci = variants
                .map { lociPosition.nucelotideLoci(it.position().toInt()) }
                .filter { it >= 0 }
                .map { it / 3 }

        val hetLociSansVariants = hetLoci subtract variantLoci

        val outputFile = "${config.outputFilePrefix}.lilac.somatic.vcf.gz"
        val vcf = LilacVCF(outputFile, config.somaticVcf)
        vcf.writeHeader()

        for (variant in variants) {
            val variantFragments = reader
                    .readFromBam(variant)
                    .map { it.qualityFilter(config.minBaseQual) }
                    .filter { it.isNotEmpty() }
                    .map { it.toAminoAcidFragment() }


            val variantFragmentAlleles = FragmentAlleles.create(variantFragments, hetLociSansVariants, winners, listOf(), listOf())
            val coverage = HlaAlleleCoverage.proteinCoverage(variantFragmentAlleles).sortedByDescending { it.totalCoverage }
            vcf.writeVariant(createNewVariant(variant, coverage))
            logger.info("    $variant -> $coverage")
        }
        vcf.close()
    }


    fun createNewVariant(variant: VariantContextDecorator, coverage: List<HlaAlleleCoverage>): VariantContext {
        val hlaLabel = if (coverage.isEmpty()) {
            listOf("UNKNOWN")
        } else {
            coverage.filter { it.totalCoverage == coverage[0].totalCoverage }.map { it.allele.toString() }
        }

        return VariantContextBuilder(variant.context()).attribute(LilacVCF.HLA, hlaLabel).make()
    }

    fun readSomatics(config: LilacConfig): List<VariantContextDecorator> {
        if (config.somaticVcf.isEmpty()) {
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