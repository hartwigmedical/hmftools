package com.hartwig.hmftools.viridian.reference;

import static java.util.Comparator.comparing;

import java.util.Comparator;

import org.jetbrains.annotations.NotNull;

// A curated virus taxon, at the rank worth telling apart for oncology, e.g. "Human papillomavirus type 16".
// It belongs to one virus group and has one or more virus genomes.
public record Virus(
        // Cleaned, human-readable name.
        String name,
        VirusGroup group
) implements Comparable<Virus>
{
    public Virus
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
    public int compareTo(Virus other)
    {
        return COMPARATOR.compare(this, other);
    }

    private static final Comparator<Virus> COMPARATOR = comparing(Virus::group).thenComparing(Virus::name);
}
