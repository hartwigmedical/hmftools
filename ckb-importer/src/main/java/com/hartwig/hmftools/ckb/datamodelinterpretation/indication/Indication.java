package com.hartwig.hmftools.ckb.datamodelinterpretation.indication;

import java.util.Date;
import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Indication {

    public abstract int id();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String source();

    @NotNull
    public abstract String definition();

    @NotNull
    public abstract String currentPreferredTerm();

    @NotNull
    public abstract Date lastUpdateDateFromDO();

    @NotNull
    public abstract List<Integer> altIds();

    @NotNull
    public abstract String termId();
}