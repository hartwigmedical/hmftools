package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Evidence {

    @Nullable
    public abstract EvidenceInfo info();

    @NotNull
    public abstract EvidenceType evidenceType();

    @Nullable
    public abstract String description();
}
