package com.hartwig.hmftools.cobalt;

import java.util.Collection;

import org.jetbrains.annotations.Nullable;

public class Chromosome
{
    public final String contig;
    public final int length;

    public Chromosome(String contig, int length)
    {
        this.contig = contig;
        this.length = length;
    }

    @Override
    public String toString() { return contig; }

    @Nullable
    public static Chromosome findByContig(String contig, Collection<Chromosome> chromosomes)
    {
        for (Chromosome c : chromosomes)
        {
            if (c.contig.equals(contig))
                return c;
        }
        return null;
    }
}
