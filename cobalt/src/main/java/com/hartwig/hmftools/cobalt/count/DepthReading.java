package com.hartwig.hmftools.cobalt.count;

import java.util.Objects;

public class DepthReading
{
    public final String Chromosome; // todo make this a Chromosome
    public final int StartPosition;
    public final double ReadDepth;
    public final double ReadGcContent;

    public DepthReading(final String chromosome, final int startPosition, final double readDepth, final double readGcContent)
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

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final DepthReading that = (DepthReading) o;
        return StartPosition == that.StartPosition && Double.compare(ReadDepth, that.ReadDepth) == 0
                && Double.compare(ReadGcContent, that.ReadGcContent) == 0 && Objects.equals(Chromosome, that.Chromosome);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(Chromosome, StartPosition, ReadDepth, ReadGcContent);
    }
}
