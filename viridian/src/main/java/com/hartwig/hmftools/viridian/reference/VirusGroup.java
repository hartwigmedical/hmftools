package com.hartwig.hmftools.viridian.reference;

import org.jetbrains.annotations.NotNull;

// A group of viruses that sets the main granularity of virus calling, e.g. "Human papillomavirus".
// Viridian reports at most one virus per group.
// It is expected that co-occurrence of virus within the same group is negligible, since Viridian can't call that.
public record VirusGroup(
        String name
) implements Comparable<VirusGroup>
{
    public VirusGroup
    {
        if(name.isEmpty())
        {
            throw new IllegalArgumentException("Invalid name");
        }
    }

    @NotNull
    @Override
    public String toString()
    {
        return name;
    }

    @Override
    public int compareTo(VirusGroup other)
    {
        return name.compareTo(other.name);
    }
}
