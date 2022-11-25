package com.hartwig.hmftools.peach.haplotype;

import com.google.common.collect.ImmutableSet;
import com.hartwig.hmftools.peach.event.HaplotypeEvent;
import org.jetbrains.annotations.NotNull;

public class WildTypeHaplotype implements Haplotype
{
    @NotNull
    public final String name;
    @NotNull
    public final ImmutableSet<HaplotypeEvent> eventsToIgnore;

    public WildTypeHaplotype(@NotNull String name, @NotNull ImmutableSet<HaplotypeEvent> eventsToIgnore)
    {
        this.name = name;
        this.eventsToIgnore = eventsToIgnore;
    }

    public boolean isRelevantFor(HaplotypeEvent event)
    {
        return false;
    }
}
