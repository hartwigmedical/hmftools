package com.hartwig.hmftools.cobalt.count;

public class ReadDepth
{
    public final String chromosome;
    public final int startPosition;
    public final double readDepth;

    public ReadDepth(final String chromosome, final int startPosition, final double readDepth)
    {
        this.chromosome = chromosome;
        this.startPosition = startPosition;
        this.readDepth = readDepth;
    }
}
