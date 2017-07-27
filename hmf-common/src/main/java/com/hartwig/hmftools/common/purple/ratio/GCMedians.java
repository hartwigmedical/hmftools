package com.hartwig.hmftools.common.purple.ratio;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

@Value.Immutable
public abstract class GCMedians implements Comparable<GCMedians> {

    public abstract int gcContent();

    public abstract int medianCount();

    @Override
    public int compareTo(@NotNull final GCMedians other) {
        return Integer.compare(gcContent(), other.gcContent());
    }
}
