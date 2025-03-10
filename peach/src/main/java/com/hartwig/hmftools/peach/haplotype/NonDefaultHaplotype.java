package com.hartwig.hmftools.peach.haplotype;

import com.google.common.collect.ImmutableList;
import com.hartwig.hmftools.peach.event.HaplotypeEvent;

public class NonDefaultHaplotype implements Haplotype
{
    private final String name;
    private final boolean isWildType;
    public final ImmutableList<HaplotypeEvent> events;

    public NonDefaultHaplotype(final String name, boolean isWildType, final ImmutableList<HaplotypeEvent> events)
    {
        if(events.isEmpty())
        {
            throw new RuntimeException(String.format("Non-wild type haplotype '%s' has no associated events", name));
        }
        this.name = name;
        this.isWildType = isWildType;
        this.events = events;
    }

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

    public boolean isRelevantFor(final HaplotypeEvent event)
    {
        return events.stream().anyMatch(event::isRelevantFor);
    }

    public int getMatchingCount(final String eventId)
    {
        return Math.toIntExact(events.stream().filter(e -> e.id().equals(eventId)).count());
    }
}
