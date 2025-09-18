package com.hartwig.hmftools.cobalt.calculations;

import java.util.Objects;

import com.hartwig.hmftools.cobalt.count.ReadDepth;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public class CobaltWindow
{
    public final Chromosome Chromosome;
    public final int Position;
    public final ReadDepth ReadDepth;
    public final GCPail GcBucket;

    public CobaltWindow(final Chromosome chromosome, final int position, final ReadDepth ReadDepth, final GCPail GcBucket)
    {
        Chromosome = chromosome;
        Position = position;
        this.ReadDepth = ReadDepth;
        this.GcBucket = GcBucket;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final CobaltWindow that = (CobaltWindow) o;
        return Position == that.Position && Chromosome == that.Chromosome;
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(Chromosome, Position);
    }
}
