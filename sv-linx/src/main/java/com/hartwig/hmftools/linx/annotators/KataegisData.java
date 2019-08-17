package com.hartwig.hmftools.linx.annotators;

public class KataegisData
{
    public final String Chromosome;
    public final long PosStart;
    public final long PosEnd;
    public final String Id;
    public final int SnvCount;

    public KataegisData(final String chromosome, long posStart, long posEnd, final String id, int snvCount)
    {
        Chromosome = chromosome;
        PosStart = posStart;
        PosEnd = posEnd;
        Id = id;
        SnvCount = snvCount;
    }

}
