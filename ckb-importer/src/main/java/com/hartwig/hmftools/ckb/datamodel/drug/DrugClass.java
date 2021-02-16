package com.hartwig.hmftools.ckb.datamodel.drug;

import java.util.Date;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DrugClass {

    public abstract int id();

    @Nullable
    public abstract Date createDate();

    @NotNull
    public abstract String drugClass();
}