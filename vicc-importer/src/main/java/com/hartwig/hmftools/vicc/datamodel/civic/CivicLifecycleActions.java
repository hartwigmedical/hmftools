package com.hartwig.hmftools.vicc.datamodel.civic;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicLifecycleActions {

    @Nullable
    public abstract CivicLastCommentedOn lastCommentedOn();

    @Nullable
    public abstract CivicLastModified lastModified();

    @Nullable
    public abstract CivicLastReviewed lastReviewed();

}
