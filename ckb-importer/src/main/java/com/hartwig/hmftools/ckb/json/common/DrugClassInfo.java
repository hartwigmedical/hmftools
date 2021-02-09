package com.hartwig.hmftools.ckb.json.common;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DrugClassInfo {

    public abstract int id();

    @NotNull
    public abstract String drugClass();
}
