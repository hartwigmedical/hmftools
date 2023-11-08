package com.hartwig.hmftools.peach.haplotype;

import com.google.common.collect.ImmutableList;
import com.hartwig.hmftools.peach.event.HaplotypeEvent;
import org.jetbrains.annotations.NotNull;

public class DefaultHaplotype implements Haplotype
{
    @NotNull
    public final String name;
    @NotNull
    public final ImmutableList<HaplotypeEvent> eventsToIgnore;
    public final boolean isWildType;

    public DefaultHaplotype(@NotNull String name, @NotNull ImmutableList<HaplotypeEvent> eventsToIgnore, boolean isWildType)
    {
        this.name = name;
        this.eventsToIgnore = eventsToIgnore;
        this.isWildType = isWildType;
    }
}
