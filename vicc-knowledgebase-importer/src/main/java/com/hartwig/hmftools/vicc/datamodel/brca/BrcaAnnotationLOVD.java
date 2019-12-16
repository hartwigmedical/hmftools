package com.hartwig.hmftools.vicc.datamodel.brca;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class BrcaAnnotationLOVD {

    @NotNull
    public abstract String variantInLOVD();

    @NotNull
    public abstract String bxId();

    @NotNull
    public abstract String dbId();

    @NotNull
    public abstract String hgvsCDNA();

    @NotNull
    public abstract String hgvsProtein();

    @NotNull
    public abstract String rna();

    @NotNull
    public abstract String variantEffect();

    @NotNull
    public abstract String variantFrequency();

    @NotNull
    public abstract String variantHaplotype();

    @NotNull
    public abstract String geneticOrigin();

    @NotNull
    public abstract String functionalAnalysisTechnique();

    @NotNull
    public abstract String functionalAnalysisResult();

    @NotNull
    public abstract String submitters();

    @NotNull
    public abstract String individuals();


}
