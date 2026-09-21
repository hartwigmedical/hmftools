package com.hartwig.hmftools.viridian.reference;

import static java.util.Comparator.comparing;

import java.util.Comparator;

import org.jetbrains.annotations.NotNull;

// TODO: maybe rename to "viral genome"? "contig" is not that descriptive. But that's a big change
// A contig which is a virus genome.
public record ViralContig(
        // Original contig name in the viral reference.
        String name,
        int length,
        // Human-readable name of the virus. Only for readability purposes.
        String virusName,
        OncologyGroup oncologyGroup
) implements Comparable<ViralContig>
{
    public ViralContig
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

    // Contig names are unique within the reference, so this orders contigs consistently with equality.
    @Override
    public int compareTo(ViralContig other)
    {
        return BY_ONCOLOGY_GROUP_THEN_NAME.compare(this, other);
    }

    private static final Comparator<ViralContig> BY_ONCOLOGY_GROUP_THEN_NAME =
            comparing(ViralContig::oncologyGroup).thenComparing(ViralContig::name);
}
