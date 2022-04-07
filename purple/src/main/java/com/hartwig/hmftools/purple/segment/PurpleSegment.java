package com.hartwig.hmftools.purple.segment;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

public class PurpleSegment implements GenomeRegion
{
    public final String Chromosome;
    public int Start;
    public int End;

    public boolean RatioSupport;
    public SegmentSupport Support;
    public boolean SvCluster;
    public int MinStart;
    public int MaxStart;

    public PurpleSegment(
            final String chromosome, final int start, final int end, final boolean ratioSupport,
            final SegmentSupport support, final boolean svCluster, final int minStart, final int maxStart)
    {
        Chromosome = chromosome;
        Start = start;
        End = end;
        RatioSupport = ratioSupport;
        Support = support;
        SvCluster = svCluster;
        MinStart = minStart;
        MaxStart = maxStart;
    }

    public static PurpleSegment from(final PurpleSegment ref)
    {
        return new PurpleSegment(
                ref.Chromosome, ref.Start, ref.End, ref.RatioSupport, ref.Support, ref.SvCluster, ref.MinStart, ref.MaxStart);
    }

    @Override
    public String chromosome()
    {
        return Chromosome;
    }

    @Override
    public int start()
    {
        return Start;
    }

    @Override
    public int end()
    {
        return End;
    }

    public String toString()
    {
        return String.format("loc(%s: %d - %d) ratio(%s) seg(%s)", Chromosome, Start, End, RatioSupport, Support);
    }
}
