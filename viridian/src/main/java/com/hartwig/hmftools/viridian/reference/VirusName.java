package com.hartwig.hmftools.viridian.reference;

import org.jetbrains.annotations.NotNull;

// The cleaned, human-readable name of a specific virus, e.g. "Human papillomavirus type 16".
// Multiple virus genomes may share one virus name, and one virus group can hold several virus names.
public record VirusName(
        String value
) implements Comparable<VirusName>
{
    public VirusName
    {
        if(value.isEmpty())
        {
            throw new IllegalArgumentException("Invalid value");
        }
    }

    @NotNull
    @Override
    public String toString()
    {
        return value;
    }

    @Override
    public int compareTo(VirusName other)
    {
        return value.compareTo(other.value);
    }
}
