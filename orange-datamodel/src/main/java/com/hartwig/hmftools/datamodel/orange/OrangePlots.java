package com.hartwig.hmftools.datamodel.orange;

import java.util.List;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface OrangePlots
{
    @Nullable
    String sageReferenceBQRPlot();

    @NotNull
    String sageTumorBQRPlot();

    @NotNull
    String purpleInputPlot();

    @NotNull
    String purpleFinalCircosPlot();

    @NotNull
    String purpleClonalityPlot();

    @NotNull
    String purpleCopyNumberPlot();

    @NotNull
    String purpleVariantCopyNumberPlot();

    @NotNull
    String purplePurityRangePlot();

    @NotNull
    String purpleKataegisPlot();

    @NotNull
    List<String> linxDriverPlots();

    @Nullable
    String cuppaSummaryPlot();

    @Nullable
    String cuppaFeaturePlot();

    @Nullable
    String cuppaChartPlot();
}
