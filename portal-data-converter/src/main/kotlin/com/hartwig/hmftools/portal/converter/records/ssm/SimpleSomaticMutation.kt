package com.hartwig.hmftools.portal.converter.records.ssm

import com.google.common.collect.Lists
import com.hartwig.hmftools.extensions.samtools.split
import com.hartwig.hmftools.portal.converter.extensions.mutationType
import com.hartwig.hmftools.portal.converter.extensions.reformatAlleles
import com.hartwig.hmftools.portal.converter.extensions.toRecord
import com.hartwig.hmftools.portal.converter.records.Record
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import org.apache.logging.log4j.LogManager
import java.io.File
import kotlin.reflect.KClass

data class SimpleSomaticMutation(private val fields: Map<SimpleSomaticMutationHeader, String>) : Record {
    companion object Factory {
        private val logger = LogManager.getLogger(SimpleSomaticMutation::class)
        val header: KClass<SimpleSomaticMutationHeader> = SimpleSomaticMutationHeader::class
        private val CHROMOSOME_STRAND = "1"                             //MIVO: 1 = forward
        private val VERIFICATION_STATUS = "2"                           //MIVO: 2 = not tested
        private operator fun invoke(analysisId: String, sampleId: String, mutationType: String, chromosome: String, start: Int, end: Int,
                                    ref: String, alt: String): SimpleSomaticMutation {
            val controlGenotype = "$ref/$ref"
            val tumorGenotype = "$alt/$alt"
            return SimpleSomaticMutation(mapOf(SimpleSomaticMutationHeader.analysis_id to analysisId,
                                               SimpleSomaticMutationHeader.analyzed_sample_id to sampleId,
                                               SimpleSomaticMutationHeader.mutation_type to mutationType,
                                               SimpleSomaticMutationHeader.chromosome to chromosome,
                                               SimpleSomaticMutationHeader.chromosome_start to start.toString(),
                                               SimpleSomaticMutationHeader.chromosome_end to end.toString(),
                                               SimpleSomaticMutationHeader.control_genotype to controlGenotype,
                                               SimpleSomaticMutationHeader.tumour_genotype to tumorGenotype,
                                               SimpleSomaticMutationHeader.reference_genome_allele to ref,
                                               SimpleSomaticMutationHeader.mutated_from_allele to ref,
                                               SimpleSomaticMutationHeader.mutated_to_allele to alt,
                                               SimpleSomaticMutationHeader.chromosome_strand to CHROMOSOME_STRAND,
                                               SimpleSomaticMutationHeader.verification_status to VERIFICATION_STATUS))
        }

        private operator fun invoke(sampleHmfId: String, variantContext: VariantContext): List<SimpleSomaticMutation> {
            val variants = if (variantContext.alternateAlleles.size > 1) {
                variantContext.split().map { it.reformatAlleles() }
            } else {
                Lists.newArrayList<VariantContext>(variantContext)
            }
            return variants.filter(this::filterVariants).map { it ->
                SimpleSomaticMutation(sampleHmfId,
                                      sampleHmfId,
                                      it.mutationType(),
                                      it.contig,
                                      it.start,
                                      it.end,
                                      it.reference.baseString,
                                      it.alternateAlleles[0].baseString)
            }
        }

        operator fun invoke(sampleHmfId: String, vcfPath: String): List<SimpleSomaticMutation> {
            val vcfReader = VCFFileReader(File(vcfPath), false)
            return vcfReader.flatMap { it -> SimpleSomaticMutation(sampleHmfId, it) }
        }

        private fun filterVariants(variant: VariantContext): Boolean {
            if (variant.mutationType() == "0") {
                logger.warn("[${variant.sampleNamesOrderedByName[0]}] could not determine mutation type for variant ${variant.contig}:${variant.start}\t${variant.reference.baseString} -> ${variant.alternateAlleles.first().baseString}")
            }
            return variant.mutationType() != "0"
        }
    }

    override fun record(): List<String> {
        return fields.toRecord()
    }
}
