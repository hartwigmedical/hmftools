package com.hartwig.hmftools.ckb.datamodel.common;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class ReferenceInfo implements ReferenceInfoExtended{

    public abstract int id();

    @Nullable
    public abstract String pubMedId();

    @Nullable
    public abstract String title();

    @Nullable
    public abstract String url();
}
