package com.hartwig.hmftools.lilac.variant

import com.hartwig.hmftools.common.drivercatalog.DriverImpact
import com.hartwig.hmftools.common.variant.VariantContextDecorator
import com.hartwig.hmftools.lilac.LilacApplication
import com.hartwig.hmftools.lilac.LilacConfig
import com.hartwig.hmftools.lilac.LociPosition
import com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage
import com.hartwig.hmftools.lilac.read.AminoAcidIndices
import com.hartwig.hmftools.lilac.read.FragmentAlleles
import com.hartwig.hmftools.lilac.read.SAMRecordReader
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci
import htsjdk.variant.vcf.VCFFileReader
import org.apache.logging.log4j.LogManager
import java.io.File

class Somatics {

    val logger = LogManager.getLogger(this::class.java)


    fun doStuff(config: LilacConfig, reader: SAMRecordReader, winners: List<HlaSequenceLoci>, lociPosition: LociPosition) {
        if (config.tumorBam.isEmpty() || config.somaticVcf.isEmpty()) {
            return
        }

        val variants = readSomatics(config)
        val variantLoci = variants
                .map { lociPosition.nucelotideLoci(it.position().toInt()) }
                .filter { it > 0 }
                .map { AminoAcidIndices.indices(it, it).first }

        for (variant in variants) {
            val alleleCount = winners.map { Pair(it.allele,0) }.toMap().toMutableMap()
            val variantFragments = reader
                    .readFromBam(variant)
                    .map { it.qualityFilter(config.minBaseQual) }
                    .filter { it.isNotEmpty() }
                    .map { it.toAminoAcidFragment() }

            val variantFragmentAlleles = variantFragments.flatMap {
                val loci = (it.aminoAcidLoci() subtract variantLoci).sorted()
                FragmentAlleles.create(listOf(it), loci, winners, listOf(), listOf())
            }

            val coverage = HlaAlleleCoverage.proteinCoverage(variantFragmentAlleles)
            logger.info("    $variant -> $coverage")
        }

    }


    fun readSomatics(config: LilacConfig): List<VariantContextDecorator> {
        if (config.somaticVcf.isEmpty()) {
            return listOf()
        }

        logger.info("Reading somatic vcf: ${config.somaticVcf}")
        val result = mutableListOf<VariantContextDecorator>()
        val fileReader = VCFFileReader(File(config.somaticVcf), false)
        for (variantContext in fileReader) {
            val enriched = VariantContextDecorator(variantContext)
            if (enriched.gene() in LilacApplication.HLA_GENES && enriched.impact() != DriverImpact.UNKNOWN && enriched.isPass) {
                result.add(enriched)
            }
        }

        logger.info("    found ${result.size} HLA variants")
        return result
    }

}