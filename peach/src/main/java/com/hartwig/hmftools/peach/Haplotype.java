package com.hartwig.hmftools.peach;

import com.google.common.collect.ImmutableSet;
import org.jetbrains.annotations.NotNull;

public class Haplotype
{
    @NotNull
    public final String gene;
    @NotNull
    public final String name;
    public final boolean wildType;
    @NotNull
    public final ImmutableSet<HaplotypeEvent> events;

    public Haplotype(@NotNull String gene, @NotNull String name, boolean wildType, @NotNull ImmutableSet<HaplotypeEvent> events)
    {
        this.gene = gene;
        this.name = name;
        this.wildType = wildType;
        this.events = events;
    }
}
