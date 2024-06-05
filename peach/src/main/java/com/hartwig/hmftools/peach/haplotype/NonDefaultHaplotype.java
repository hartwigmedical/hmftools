package com.hartwig.hmftools.peach.haplotype;

import com.google.common.collect.ImmutableList;
import com.hartwig.hmftools.peach.event.HaplotypeEvent;

import org.jetbrains.annotations.NotNull;

public class NonDefaultHaplotype implements Haplotype
{
    @NotNull
    private final String name;
    private final boolean isWildType;
    @NotNull
    public final ImmutableList<HaplotypeEvent> events;

    public NonDefaultHaplotype(@NotNull String name, boolean isWildType, @NotNull ImmutableList<HaplotypeEvent> events)
    {
        if(events.isEmpty())
        {
            throw new RuntimeException(String.format("Non-wild type haplotype '%s' has no associated events", name));
        }
        this.name = name;
        this.isWildType = isWildType;
        this.events = events;
    }

    @NotNull
    @Override
    public String getName()
    {
        return name;
    }

    @Override
    public boolean isWildType()
    {
        return isWildType;
    }

    public boolean isRelevantFor(@NotNull HaplotypeEvent event)
    {
        return events.stream().anyMatch(event::isRelevantFor);
    }

    public int getMatchingCount(@NotNull String eventId)
    {
        return Math.toIntExact(events.stream().filter(e -> e.id().equals(eventId)).count());
    }
}
