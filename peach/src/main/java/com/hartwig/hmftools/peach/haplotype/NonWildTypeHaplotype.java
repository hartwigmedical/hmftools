package com.hartwig.hmftools.peach.haplotype;

import com.google.common.collect.ImmutableList;
import com.hartwig.hmftools.peach.event.HaplotypeEvent;
import org.jetbrains.annotations.NotNull;

public class NonWildTypeHaplotype implements Haplotype
{
    @NotNull
    public final String name;
    @NotNull
    public final ImmutableList<HaplotypeEvent> events;

    public NonWildTypeHaplotype(@NotNull String name, @NotNull ImmutableList<HaplotypeEvent> events)
    {
        if (events.size() < 1)
            throw new RuntimeException(String.format("Non-wild type haplotype '%s' has no associated events", name));
        this.name = name;
        this.events = events;
    }

    public boolean isRelevantFor(HaplotypeEvent event)
    {
        return events.stream().anyMatch(event::isRelevantFor);
    }
}
