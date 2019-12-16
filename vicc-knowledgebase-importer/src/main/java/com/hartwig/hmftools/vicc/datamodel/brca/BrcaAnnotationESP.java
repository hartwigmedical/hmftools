package com.hartwig.hmftools.vicc.datamodel.brca;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class BrcaAnnotationESP {

    @NotNull
    public abstract String variantInESP();

    @NotNull
    public abstract String bxId();

    @NotNull
    public abstract String minorAlleleFrequencyPercent();

    @NotNull
    public abstract String alleleFrequency();

    @NotNull
    public abstract String aaAlleleFrequency();

    @NotNull
    public abstract String eaAlleleFrequency();

}
