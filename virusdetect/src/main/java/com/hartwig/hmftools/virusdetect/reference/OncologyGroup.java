package com.hartwig.hmftools.virusdetect.reference;

import org.jetbrains.annotations.NotNull;

// A group of viruses at the level of taxonomy granularity which matters to us.
public record OncologyGroup(String name)
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
}
