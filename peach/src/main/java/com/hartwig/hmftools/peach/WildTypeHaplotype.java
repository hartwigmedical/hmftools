package com.hartwig.hmftools.peach;

import com.google.common.collect.ImmutableSet;
import org.jetbrains.annotations.NotNull;

public class WildTypeHaplotype implements Haplotype
{
    @NotNull
    public final String name;
    @NotNull
    public final ImmutableSet<HaplotypeEvent> ignoredEvents;

    public WildTypeHaplotype(@NotNull String name, @NotNull ImmutableSet<HaplotypeEvent> ignoredEvents)
    {
        this.name = name;
        this.ignoredEvents = ignoredEvents;
    }

    public boolean isRelevantFor(HaplotypeEvent event)
    {
        return false;
    }
}
