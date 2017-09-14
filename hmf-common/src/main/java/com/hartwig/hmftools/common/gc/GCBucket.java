package com.hartwig.hmftools.common.gc;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(of = "new",
             allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public interface GCBucket {

    @NotNull
    static GCBucket create(@NotNull final GCProfile profile) {
        return new ImmutableGCBucket((int) Math.round(profile.gcContent() * 100));
    }

    int bucket();
}
