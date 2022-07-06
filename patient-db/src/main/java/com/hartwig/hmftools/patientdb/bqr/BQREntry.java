package com.hartwig.hmftools.patientdb.bqr;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class BQREntry {

    @NotNull
    public abstract String ref();

    @NotNull
    public abstract String alt();

    @NotNull
    public abstract String trinucleotideContext();

    public abstract int count();

    public abstract double origQuality();

    public abstract double recalibratedQuality();

}
