package com.hartwig.hmftools.ckb.datamodelinterpretation.therapy;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Reference {

    public abstract int id();

    @NotNull
    public abstract String pudMedId();

    @NotNull
    public abstract String title();

    @NotNull
    public abstract String url();
}
