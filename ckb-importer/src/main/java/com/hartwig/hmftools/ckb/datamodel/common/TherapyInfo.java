package com.hartwig.hmftools.ckb.datamodel.common;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class TherapyInfo {

    public abstract int id();

    @NotNull
    public abstract String therapyName();

    @Nullable
    public abstract String synonyms();
}
