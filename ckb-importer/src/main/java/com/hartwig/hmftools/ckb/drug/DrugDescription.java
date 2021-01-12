package com.hartwig.hmftools.ckb.drug;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DrugDescription {

    @NotNull
    public abstract String description();

    @NotNull
    public abstract List<DrugReference> references();
}
