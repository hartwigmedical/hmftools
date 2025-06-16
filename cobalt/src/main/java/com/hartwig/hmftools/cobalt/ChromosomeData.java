package com.hartwig.hmftools.cobalt;

public class ChromosomeData
{
    public final String Contig;
    public final int Length;

    public ChromosomeData(String contig, int length)
    {
        Contig = contig;
        Length = length;
    }

    @Override
    public String toString() { return Contig; }
}
