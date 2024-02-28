package com.hartwig.hmftools.peach.haplotype;

import com.google.common.collect.ImmutableList;
import com.hartwig.hmftools.peach.event.HaplotypeEvent;

import org.jetbrains.annotations.NotNull;

public class DefaultHaplotype implements Haplotype
{
    @NotNull
    private final String name;
    private final boolean isWildType;
    @NotNull
    public final ImmutableList<HaplotypeEvent> eventsToIgnore;

    public DefaultHaplotype(@NotNull String name, boolean isWildType, @NotNull ImmutableList<HaplotypeEvent> eventsToIgnore)
    {
        this.name = name;
        this.isWildType = isWildType;
        this.eventsToIgnore = eventsToIgnore;
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
}
