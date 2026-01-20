package com.hartwig.hmftools.purple.segment;

import java.util.List;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.purple.SegmentSupport;

public class PurpleSupportSegment implements GenomeRegion
{
    public final String Chromosome;
    public final boolean RatioSupport;

    public SegmentSupport Support;
    public boolean SvCluster;

    private int mStart;
    private int mEnd;
    private int mMinStart;
    private int mMaxStart;

    public PurpleSupportSegment(
            final String chromosome, final int start, final int end, final boolean ratioSupport,
            final SegmentSupport support, final boolean svCluster)
    {
        this(chromosome, start, end, ratioSupport, support, svCluster, start, start);
    }

    public PurpleSupportSegment(
            final String chromosome, final int start, final int end, final boolean ratioSupport,
            final SegmentSupport support, final boolean svCluster, final int minStart, final int maxStart)
    {
        Chromosome = chromosome;
        mStart = start;
        mEnd = end;
        RatioSupport = ratioSupport;
        Support = support;
        SvCluster = svCluster;
        mMinStart = minStart;
        mMaxStart = maxStart;
    }

    public static PurpleSupportSegment from(final PurpleSupportSegment ref)
    {
        return new PurpleSupportSegment(
                ref.Chromosome, ref.mStart, ref.mEnd, ref.RatioSupport, ref.Support, ref.SvCluster, ref.mMinStart, ref.mMaxStart);
    }

    @Override
    public String chromosome()
    {
        return Chromosome;
    }

    @Override
    public int start()
    {
        return mStart;
    }

    @Override
    public int end()
    {
        return mEnd;
    }

    public void setStart(int position)
    {
        mStart = position;
    }

    public void setEnd(int position)
    {
        mEnd = position;
    }

    public void setMinStart(int position)
    {
        mMinStart = position;
    }

    public void setMaxStart(int position)
    {
        mMaxStart = position;
    }

    public int minStart()
    {
        return mMinStart;
    }

    public int maxStart()
    {
        return mMaxStart;
    }

    public List<PurpleSupportSegment> split(int position)
    {
        if(position <= mStart || position > mEnd)
        {
            return List.of(this);
        }
        int endOfLeft = position - 1;
        int newMaxStart = Math.min(mMaxStart, endOfLeft);
        return List.of(new PurpleSupportSegment(Chromosome, mStart, endOfLeft, RatioSupport, Support, SvCluster, mMinStart, newMaxStart),
                new PurpleSupportSegment(Chromosome, position, mEnd, RatioSupport, Support, SvCluster));
    }

    public List<PurpleSupportSegment> splitBy(GenomeRegion other)
    {
        if(!chromosome().equals(other.chromosome()))
        {
            return List.of(this);
        }

        return split(other.start()).stream()
                .flatMap(segment -> segment.split(other.end() + 1).stream())
                .toList();
    }

    public String toString()
    {
        return String.format("loc(%s: %d - %d) ratio(%s) seg(%s)", Chromosome, mStart, mEnd, RatioSupport, Support);
    }
}
