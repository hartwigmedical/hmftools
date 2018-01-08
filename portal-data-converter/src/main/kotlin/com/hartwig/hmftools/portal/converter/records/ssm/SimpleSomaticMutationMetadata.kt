package com.hartwig.hmftools.portal.converter.records.ssm

import com.hartwig.hmftools.portal.converter.extensions.toRecord
import com.hartwig.hmftools.portal.converter.records.Record
import kotlin.reflect.KClass

data class SimpleSomaticMutationMetadata(private val fields: Map<SimpleSomaticMutationMetadataHeader, String>) : Record {
    companion object Factory {
        val header: KClass<SimpleSomaticMutationMetadataHeader> = SimpleSomaticMutationMetadataHeader::class
        private val ASSEMBLY_VERSION = "1"                 // GRCh37
        private val PLATFORM = "84"                        // Illumina HiSeq X Ten
        private val SEQUENCING_STRATEGY = "1"              // WGS
        private val ALIGNMENT_ALGORITHM = "BWA-MEM"
        private val CALLING_ALGORITHM = "Strelka"
        private val RAW_DATA_ACCESSION = "HMF0"
        private val RAW_DATA_REPOSITORY = "1"              // EGA; TODO: replace
        operator fun invoke(sampleId: String): SimpleSomaticMutationMetadata {
            return SimpleSomaticMutationMetadata(mapOf(
                    SimpleSomaticMutationMetadataHeader.analysis_id to sampleId,
                    SimpleSomaticMutationMetadataHeader.analyzed_sample_id to sampleId,
                    SimpleSomaticMutationMetadataHeader.assembly_version to ASSEMBLY_VERSION,
                    SimpleSomaticMutationMetadataHeader.platform to PLATFORM,
                    SimpleSomaticMutationMetadataHeader.sequencing_strategy to SEQUENCING_STRATEGY,
                    SimpleSomaticMutationMetadataHeader.alignment_algorithm to ALIGNMENT_ALGORITHM,
                    SimpleSomaticMutationMetadataHeader.variation_calling_algorithm to CALLING_ALGORITHM,
                    SimpleSomaticMutationMetadataHeader.raw_data_accession to RAW_DATA_ACCESSION,
                    SimpleSomaticMutationMetadataHeader.raw_data_repository to RAW_DATA_REPOSITORY
            ))
        }
    }

    override fun record(): List<String> {
        return fields.toRecord()
    }
}