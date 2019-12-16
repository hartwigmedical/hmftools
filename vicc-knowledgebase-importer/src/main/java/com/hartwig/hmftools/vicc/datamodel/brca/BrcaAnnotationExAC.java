package com.hartwig.hmftools.vicc.datamodel.brca;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class BrcaAnnotationExAC {

    @NotNull
    public abstract String variantInExAC();

    @NotNull
    public abstract String bxId();

    @NotNull
    public abstract String alleleFrequency();

    @NotNull
    public abstract String alleleFrequencyAFR();

    @NotNull
    public abstract String alleleFrequencyAMR();

    @NotNull
    public abstract String alleleFrequencyEAS();

    @NotNull
    public abstract String alleleFrequencyFIN();

    @NotNull
    public abstract String alleleFrequencyNFE();

    @NotNull
    public abstract String alleleFrequencyOTH();

    @NotNull
    public abstract String alleleFrequencySAS();

    @NotNull
    public abstract String alleleNumberAFR();

    @NotNull
    public abstract String alleleNumberAMR();

    @NotNull
    public abstract String alleleNumberEAS();

    @NotNull
    public abstract String alleleNumberFIN();

    @NotNull
    public abstract String alleleNumberNFE();

    @NotNull
    public abstract String alleleNumberOTH();

    @NotNull
    public abstract String alleleNumberSAS();

    @NotNull
    public abstract String homozygousCountAFR();

    @NotNull
    public abstract String homozygousCountAMR();

    @NotNull
    public abstract String homozygousCountEAS();

    @NotNull
    public abstract String homozygousCountFIN();

    @NotNull
    public abstract String homozygousCountNFE();

    @NotNull
    public abstract String homozygousCountOTH();

    @NotNull
    public abstract String homozygousCountSAS();

    @NotNull
    public abstract String alleleCountAFR();

    @NotNull
    public abstract String alleleCountAMR();

    @NotNull
    public abstract String alleleCountEAS();

    @NotNull
    public abstract String alleleCountFIN();

    @NotNull
    public abstract String alleleCountNFE();

    @NotNull
    public abstract String alleleCountOTH();

    @NotNull
    public abstract String alleleCountSAS();




}
