package com.hartwig.hmftools.vicc.datamodel.brca;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class BrcaAnnotation1000Genomes {

    @NotNull
    public abstract String variantIn1000Genomes();

    @NotNull
    public abstract String bxId();

    @NotNull
    public abstract String alleleFrequency();

    @NotNull
    public abstract String afrAlleleFrequency();

    @NotNull
    public abstract String amrAlleleFrequency();

    @NotNull
    public abstract String easAlleleFrequency();

    @NotNull
    public abstract String eurAlleleFrequency();

    @NotNull
    public abstract String sasAlleleFrequency();
}
