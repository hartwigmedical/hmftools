package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicLifecycleActions {

    @NotNull
    public abstract CivicLastCommentedOn lastCommentedOn();

    @NotNull
    public abstract CivicLastModified lastModified();

    @NotNull
    public abstract CivicLastReviewed lastReviewed();



}
