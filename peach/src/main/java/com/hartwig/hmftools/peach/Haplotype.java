package com.hartwig.hmftools.peach;

import com.google.common.collect.ImmutableSet;
import org.jetbrains.annotations.NotNull;

import java.util.Set;

public class Haplotype
{
    @NotNull
    public final String name;
    public final boolean wildType;
    @NotNull
    public final ImmutableSet<HaplotypeEvent> events;

    public Haplotype(@NotNull String name, boolean wildType, @NotNull ImmutableSet<HaplotypeEvent> events)
    {
        this.name = name;
        this.wildType = wildType;
        this.events = events;
    }

    public boolean isRelevantFor(HaplotypeEvent event)
    {
        if (event instanceof VariantHaplotypeEvent)
        {
            return isRelevantFor((VariantHaplotypeEvent) event);
        }
        else
        {
            throw new RuntimeException(String.format("Cannot determine relevance of event %s", event));
        }
    }

    private boolean isRelevantFor(VariantHaplotypeEvent event)
    {
        return events.stream()
                .filter(e -> e instanceof VariantHaplotypeEvent)
                .map(e -> (VariantHaplotypeEvent) e)
                .map(VariantHaplotypeEvent::getCoveredPositions)
                .flatMap(Set::stream)
                .anyMatch(event.getCoveredPositions()::contains);
    }
}
