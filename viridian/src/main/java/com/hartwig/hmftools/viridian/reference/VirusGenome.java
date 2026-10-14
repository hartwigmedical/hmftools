package com.hartwig.hmftools.viridian.reference;

import static java.util.Comparator.comparing;

import java.util.Comparator;

import org.jetbrains.annotations.NotNull;

// A contig which is a virus genome.
public record VirusGenome(
        // Original contig name in the virus reference.
        String name,
        int length,
        // Human-readable name of the virus. Only for readability purposes.
        String virusName,
        OncologyGroup oncologyGroup
) implements Comparable<VirusGenome>
{
    public VirusGenome
    {
        if(length < 1)
        {
            throw new IllegalArgumentException("Invalid contig length: " + length);
        }
    }

    @NotNull
    @Override
    public String toString()
    {
        return name;
    }

    @Override
    public int compareTo(@NotNull VirusGenome other)
    {
        Comparator<VirusGenome> comparator =
                comparing(VirusGenome::oncologyGroup).thenComparing(VirusGenome::name);
        return comparator.compare(this, other);
    }
}
