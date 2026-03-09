package com.hartwig.hmftools.finding.datamodel;

import java.util.Set;

import org.jspecify.annotations.Nullable;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record HlaAllele(
        @NotNull String findingKey,
        @NotNull String geneClass,
        @NotNull String gene,
        @NotNull String allele,
        @NotNull String alleleGroup,
        @NotNull String hlaProtein,
        @NotNull Set<QcStatus> qcStatus,
        int germlineCopyNumber,
        @Nullable Double tumorCopyNumber,
        @Nullable Integer refFragments,
        int tumorFragments,
        @Nullable Integer rnaFragments,
        double somaticMissense,
        double somaticNonsenseOrFrameshift,
        double somaticSplice,
        double somaticSynonymous,
        double somaticInframeIndel
) implements Finding
{
    public enum QcStatus
    {
        PASS,
        FAIL,
        WARN_UNMATCHED_ALLELE,
        WARN_UNMATCHED_SOMATIC_VARIANT,
        WARN_UNMATCHED_HAPLOTYPE,
        WARN_UNMATCHED_AMINO_ACID,
        WARN_LOW_COVERAGE,
        WARN_LOW_BASE_QUAL,
        WARN_UNMATCHED_INDEL
    }

    public boolean hasSomaticVariants()
    {
        return Doubles.positive(somaticMissense()) || Doubles.positive(somaticNonsenseOrFrameshift()) || Doubles.positive(
                somaticSplice()) || Doubles.positive(somaticInframeIndel());
    }
}
