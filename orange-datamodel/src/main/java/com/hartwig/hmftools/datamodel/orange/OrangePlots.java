package com.hartwig.hmftools.datamodel.orange;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.List;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class OrangePlots {

    @Nullable
    public abstract String sageReferenceBQRPlot();

    @NotNull
    public abstract String sageTumorBQRPlot();

    @NotNull
    public abstract String purpleInputPlot();

    @NotNull
    public abstract String purpleFinalCircosPlot();

    @NotNull
    public abstract String purpleClonalityPlot();

    @NotNull
    public abstract String purpleCopyNumberPlot();

    @NotNull
    public abstract String purpleVariantCopyNumberPlot();

    @NotNull
    public abstract String purplePurityRangePlot();

    @Nullable
    public abstract String purpleKataegisPlot();

    @NotNull
    public abstract List<String> linxDriverPlots();

    @Nullable
    public abstract String cuppaSummaryPlot();

    @Nullable
    public abstract String cuppaFeaturePlot();

    @Nullable
    public abstract String cuppaChartPlot();
}
