package com.hartwig.hmftools.viridian.integration.variant_extract;

import org.jetbrains.annotations.Nullable;

// ESVEE's read support at one breakend.
public record BreakendSupport(
        int tumorVariantFragments,
        int tumorReferenceFragments,
        double tumorAlleleFrequency,
        // Absent for a tumor-only VCF.
        @Nullable Integer normalVariantFragments,
        @Nullable Integer normalReferenceFragments
)
{
}
