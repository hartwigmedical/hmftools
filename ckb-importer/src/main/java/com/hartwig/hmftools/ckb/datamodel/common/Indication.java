package com.hartwig.hmftools.ckb.datamodel.common;

import java.util.Date;
import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Indication {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String source();

    @Nullable
    public abstract String definition();

    @Nullable
    public abstract String currentPreferredTerm();

    @Nullable
    public abstract Date lastUpdateDateFromDO();

    @NotNull
    public abstract List<String> altIds();

    @NotNull
    public abstract String termId();
}