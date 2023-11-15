package com.hartwig.hmftools.cobalt.count;

public class ReadDepth
{
    public final String chromosome;
    public final int startPosition;
    public final double readDepth;
    public final double readGcContent;

    public ReadDepth(final String chromosome, final int startPosition, final double readDepth, final double readGcContent)
    {
        this.chromosome = chromosome;
        this.startPosition = startPosition;
        this.readDepth = readDepth;
        this.readGcContent = readGcContent;
    }
}
