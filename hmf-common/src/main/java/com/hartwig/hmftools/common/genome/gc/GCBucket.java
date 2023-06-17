package com.hartwig.hmftools.common.genome.gc;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(of = "new",
             allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public interface GCBucket
{
    static GCBucket create(final GCProfile profile)
    {
        return new ImmutableGCBucket(calcGcBucket(profile.gcContent() * 100));
    }

    int bucket();

    static int calcGcBucket(double gcContent) { return (int) Math.round(gcContent * 100); }
}
