package com.hartwig.hmftools.peach.haplotype;

import com.google.common.collect.ImmutableList;
import com.hartwig.hmftools.peach.event.HaplotypeEvent;
import org.jetbrains.annotations.NotNull;

public class NonDefaultHaplotype implements Haplotype
{
    @NotNull
    public final String name;
    @NotNull
    public final ImmutableList<HaplotypeEvent> events;
    public final boolean isWildType;

    public NonDefaultHaplotype(@NotNull String name, @NotNull ImmutableList<HaplotypeEvent> events, boolean isWildType)
    {
        if (events.isEmpty())
            throw new RuntimeException(String.format("Non-wild type haplotype '%s' has no associated events", name));
        this.name = name;
        this.events = events;
        this.isWildType = isWildType;
    }

    public boolean isRelevantFor(HaplotypeEvent event)
    {
        return events.stream().anyMatch(event::isRelevantFor);
    }
}
