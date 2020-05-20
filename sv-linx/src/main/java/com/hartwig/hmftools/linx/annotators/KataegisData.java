package com.hartwig.hmftools.linx.annotators;

public class KataegisData
{
    public final String Chromosome;
    public final int PosStart;
    public final int PosEnd;
    public final String Id;
    public final int SnvCount;

    public KataegisData(final String chromosome, int posStart, int posEnd, final String id, int snvCount)
    {
        Chromosome = chromosome;
        PosStart = posStart;
        PosEnd = posEnd;
        Id = id;
        SnvCount = snvCount;
    }

}
