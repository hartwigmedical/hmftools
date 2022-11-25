package com.hartwig.hmftools.peach.haplotype;

import com.google.common.collect.ImmutableList;
import com.hartwig.hmftools.peach.event.HaplotypeEvent;
import com.hartwig.hmftools.peach.event.VariantHaplotypeEvent;
import org.jetbrains.annotations.NotNull;

import java.util.Set;

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
