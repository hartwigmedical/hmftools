package com.hartwig.hmftools.finding.datamodel;

import java.util.List;

import org.jspecify.annotations.Nullable;

import io.soabase.recordbuilder.core.RecordBuilder;
import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record VisualisationFiles(
        @NotNull String purpleInputPlot,
        @NotNull String purpleFinalCircosPlot,
        @NotNull String purpleClonalityPlot,
        @NotNull String purpleCopyNumberPlot,
        @NotNull String purpleVariantCopyNumberPlot,
        @NotNull String purplePurityRangePlot,
        @NotNull String purpleKataegisPlot,
        @NotNull String qseePlot,
        @NotNull List<String> linxDriverPlots,
        @NotNull List<String> sageVisualisations,
        @Nullable String cuppaSummaryPlot)
{
}
