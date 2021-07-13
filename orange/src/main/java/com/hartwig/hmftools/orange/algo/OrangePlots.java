package com.hartwig.hmftools.orange.algo;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class OrangePlots {

    @NotNull
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

    @NotNull
    public abstract String purpleKataegisPlot();

    @NotNull
    public abstract List<String> linxDriverPlots();
}
