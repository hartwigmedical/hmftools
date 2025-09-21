package com.hartwig.hmftools.cobalt.count;

public class ReadDepth
{
    public final String Chromosome; // todo make this a Chromosome
    public final int StartPosition;
    public final double ReadDepth;
    public final double ReadGcContent;

    public ReadDepth(final String chromosome, final int startPosition, final double readDepth, final double readGcContent)
    {
        Chromosome = chromosome;
        StartPosition = startPosition;
        ReadDepth = readDepth;
        ReadGcContent = readGcContent;
    }

    @Override
    public String toString()
    {
        return "ReadDepth{" +
                "Chromosome='" + Chromosome + '\'' +
                ", StartPosition=" + StartPosition +
                ", ReadDepth=" + ReadDepth +
                ", ReadGcContent=" + ReadGcContent +
                '}';
    }
}
