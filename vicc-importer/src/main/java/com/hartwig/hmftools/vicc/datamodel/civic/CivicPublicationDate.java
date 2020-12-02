package com.hartwig.hmftools.vicc.datamodel.civic;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicPublicationDate {

    @Nullable
    public abstract String year();

    @Nullable
    public abstract String month();

    @Nullable
    public abstract String day();
}
