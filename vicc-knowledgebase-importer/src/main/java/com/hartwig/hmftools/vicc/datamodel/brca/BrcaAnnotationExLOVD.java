package com.hartwig.hmftools.vicc.datamodel.brca;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class BrcaAnnotationExLOVD {

    @NotNull
    public abstract String variantInExLOVD();

    @NotNull
    public abstract String bxId();

    @NotNull
    public abstract String cooccurrenceLR();

    @NotNull
    public abstract String sumFamilyLR();

    @NotNull
    public abstract String segregationLR();

    @NotNull
    public abstract String posteriorProbability();

    @NotNull
    public abstract String missenseAnalysisPriorProbability();

    @NotNull
    public abstract String combinedPriorProbability();

    @NotNull
    public abstract String iarcClass();

    @NotNull
    public abstract String literatureSource();
}
