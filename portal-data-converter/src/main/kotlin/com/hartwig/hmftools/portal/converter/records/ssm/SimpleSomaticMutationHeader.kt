package com.hartwig.hmftools.portal.converter.records.ssm

enum class SimpleSomaticMutationHeader {
    analysis_id,
    analyzed_sample_id,
    mutation_type,
    chromosome,
    chromosome_start,
    chromosome_end,
    chromosome_strand,
    reference_genome_allele,
    control_genotype,
    mutated_from_allele,
    mutated_to_allele,
    tumour_genotype,
    expressed_allele,
    quality_score,
    probability,
    total_read_count,
    mutant_allele_read_count,
    verification_status,
    verification_platform,
    biological_validation_status,
    biological_validation_platform
}
