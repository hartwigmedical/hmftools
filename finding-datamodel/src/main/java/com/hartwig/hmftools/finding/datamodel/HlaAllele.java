package com.hartwig.hmftools.finding.datamodel;

import org.jspecify.annotations.Nullable;

import io.soabase.recordbuilder.core.RecordBuilder;
import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record HlaAllele(
        @NotNull String findingKey,
        @NotNull String event,
        @NotNull String gene,
        @NotNull String allele,
        @NotNull String alleleGroup,
        @NotNull String hlaProtein,
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
) implements Event
{
    public boolean hasSomaticVariants()
    {
        return Doubles.positive(somaticMissense()) || Doubles.positive(somaticNonsenseOrFrameshift()) || Doubles.positive(
                somaticSplice()) || Doubles.positive(somaticInframeIndel());
    }
}
