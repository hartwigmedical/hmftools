package com.hartwig.hmftools.finding.datamodel;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record PurityPloidyFit(
        @NotNull FittedPurityMethod fittedPurityMethod,
        double purity,
        double minPurity,
        double maxPurity,
        double ploidy,
        double minPloidy,
        double maxPloidy
)
{
    public enum FittedPurityMethod
    {
        NORMAL,
        HIGHLY_DIPLOID,
        SOMATIC,
        NO_TUMOR
    }
}
