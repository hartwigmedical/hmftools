package com.hartwig.hmftools.ckb.common;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DescriptionInfo {

    @Nullable
    public abstract String description();

    @NotNull
    public abstract List<ReferenceInfo> reference();
}
