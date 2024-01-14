package com.hartwig.hmftools.cobalt;

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
}
