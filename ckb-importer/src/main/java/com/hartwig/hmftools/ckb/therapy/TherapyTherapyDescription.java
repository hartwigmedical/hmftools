package com.hartwig.hmftools.ckb.therapy;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class TherapyTherapyDescription {

    @Nullable
    public abstract String description();

    @NotNull
    public abstract List<TherapyReference> reference();


}
