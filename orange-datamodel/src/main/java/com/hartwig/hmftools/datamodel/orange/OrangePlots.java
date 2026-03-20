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
    @NotNull
    String purpleInputCircosPlot();

    @NotNull
    String purpleFinalCircosPlot();

    @NotNull
    String purpleClonalityPlot();

    @NotNull
    String purpleCopyNumberPlot();

    @NotNull
    String purpleMinorAlleleMapPlot();

    @NotNull
    String purpleVariantCopyNumberPlot();

    @NotNull
    String purplePurityRangePlot();

    @NotNull
    String purpleRainfallPlot();

    @NotNull
    List<String> linxDriverPlots();

    @Nullable
    String cuppaSummaryPlot();

    @Nullable
    String qSeePlot();
}
