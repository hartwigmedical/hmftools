package com.hartwig.hmftools.viridian.reference;

import static java.util.Comparator.comparing;

import java.util.Comparator;

import org.jetbrains.annotations.NotNull;

// A virus genome in the virus reference, with a FASTA contig.
public record VirusGenome(
        // Contig name in the FASTA and matching the original virus database.
        String contigName,
        int length,
        Virus virus
) implements Comparable<VirusGenome>
{
    public VirusGenome
    {
        if(length < 1)
        {
            throw new IllegalArgumentException("Invalid virus genome length: " + length);
        }
    }

    public VirusGroup virusGroup()
    {
        return virus.group();
    }

    @NotNull
    @Override
    public String toString()
    {
        return contigName;
    }

    @Override
    public int compareTo(VirusGenome other)
    {
        return COMPARATOR.compare(this, other);
    }

    private static final Comparator<VirusGenome> COMPARATOR =
            comparing(VirusGenome::virusGroup).thenComparing(VirusGenome::contigName);
}
