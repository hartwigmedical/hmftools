package com.hartwig.hmftools.virusdetect.integration;

import org.jetbrains.annotations.Nullable;

// ESVEE's read support at one breakend. Each field is null when the VCF carries no genotype for that sample.
public record BreakendSupport(
        @Nullable Integer tumorVariantFragments,
        @Nullable Integer tumorReferenceFragments,
        @Nullable Double tumorAlleleFrequency,
        @Nullable Integer normalVariantFragments,
        @Nullable Integer normalReferenceFragments
)
{
}
