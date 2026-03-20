package com.hartwig.hmftools.finding.datamodel;

import jakarta.validation.constraints.NotNull;

@SuppressWarnings("unused")
@RecordBuilder
public record PurityPloidyFit(
        @NotNull FittedPurityMethod fittedPurityMethod,
        double purity,
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
        @NotNull VisualisationFile purpleRainfallPlot
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
