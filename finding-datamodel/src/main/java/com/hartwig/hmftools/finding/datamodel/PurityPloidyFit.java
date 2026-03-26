package com.hartwig.hmftools.finding.datamodel;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record PurityPloidyFit(
        @NotNull FittedPurityMethod fittedPurityMethod,
        @NotNull ThresholdValue purity,
        double minPurity,
        double maxPurity,
        double ploidy,
        double minPloidy,
        double maxPloidy,
        @NotNull VisualisationFile purpleInputPlot,
        @NotNull VisualisationFile purpleCircosPlot,
        @NotNull VisualisationFile purpleClonalityPlot,
        @NotNull VisualisationFile purpleCopyNumberPlot,
        @NotNull VisualisationFile purpleVariantCopyNumberPlot,
        @NotNull VisualisationFile purplePurityRangePlot,
        @NotNull VisualisationFile purpleKataegisPlot
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
