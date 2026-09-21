package com.hartwig.hmftools.viridian.reference;

import org.jetbrains.annotations.NotNull;

// A group of viruses at the level of taxonomy granularity which matters to us.
public record OncologyGroup(
    String name
) implements Comparable<OncologyGroup>
{
    public OncologyGroup
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
    public int compareTo(OncologyGroup other)
    {
        return name.compareTo(other.name);
    }
}
