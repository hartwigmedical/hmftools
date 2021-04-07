package com.hartwig.hmftools.lilac.variant

import com.hartwig.hmftools.common.drivercatalog.DriverImpact
import com.hartwig.hmftools.common.variant.VariantContextDecorator
import com.hartwig.hmftools.lilac.LilacApplication
import com.hartwig.hmftools.lilac.LilacConfig
import com.hartwig.hmftools.lilac.LociPosition
import com.hartwig.hmftools.lilac.read.AminoAcidIndices
import com.hartwig.hmftools.lilac.read.SAMRecordReader
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci
import htsjdk.variant.vcf.VCFFileReader
import org.apache.logging.log4j.LogManager
import java.io.File

class Somatics {

    val logger = LogManager.getLogger(this::class.java)


    fun doStuff(config: LilacConfig, reader: SAMRecordReader, winners: List<HlaSequenceLoci>, lociPosition: LociPosition) {
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

            for (variantFragment in variantFragments) {
                val loci = (variantFragment.aminoAcidLoci() subtract variantLoci).sorted().toIntArray()
                val lociSequence = variantFragment.aminoAcids(*loci)

                for (winner in winners) {
                    if (winner.consistentWith(lociSequence, *loci)) {
                        alleleCount[winner.allele] = alleleCount[winner.allele]!! + 1
                    }
                }
            }

            println(variant)
            println(alleleCount)
            println("Sdf")
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