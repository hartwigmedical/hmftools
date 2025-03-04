package com.hartwig.hmftools.peach.haplotype;

import com.google.common.collect.ImmutableList;
import com.hartwig.hmftools.peach.event.HaplotypeEvent;

public class DefaultHaplotype implements Haplotype
{
    private final String name;
    private final boolean isWildType;
    public final ImmutableList<HaplotypeEvent> eventsToIgnore;

    public DefaultHaplotype(final String name, boolean isWildType, final ImmutableList<HaplotypeEvent> eventsToIgnore)
    {
        this.name = name;
        this.isWildType = isWildType;
        this.eventsToIgnore = eventsToIgnore;
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
}
