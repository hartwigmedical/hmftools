package com.hartwig.hmftools.finding.datamodel;

import java.util.List;
import java.util.Map;

import org.jspecify.annotations.Nullable;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record VisualisationFiles(
        @Nullable String referenceBqrPlot,
        @NotNull String tumorBqrPlot,
        @NotNull String purpleInputPlot,
        @NotNull String purpleFinalCircosPlot,
        @NotNull String purpleClonalityPlot,
        @NotNull String purpleCopyNumberPlot,
        @NotNull String purpleVariantCopyNumberPlot,
        @NotNull String purplePurityRangePlot,
        @NotNull String purpleKataegisPlot,
        @NotNull String qseePlot,
        @NotNull List<String> linxDriverPlots,
        @NotNull Map<String, String> sageVisualisations,
        @Nullable String cuppaSummaryPlot)
{
}
