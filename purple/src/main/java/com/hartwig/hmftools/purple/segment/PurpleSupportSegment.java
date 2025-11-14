package com.hartwig.hmftools.purple.segment;

import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

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

    public void setStart(int position)  { mStart = position; }

    public void setEnd(int position)
    {
        mEnd = position;

        // if(checkPositions)
        //    checkPositions();
    }

    public void setMinStart(int position) { mMinStart = position; }
    public void setMaxStart(int position) { mMaxStart = position; }
    public int minStart()
    {
        return mMinStart;
    }
    public int maxStart() { return mMaxStart; }

    private void checkPositions()
    {
        if(mEnd < mStart || mMaxStart > mEnd || mMaxStart > mEnd)
        {
            PPL_LOGGER.warn("invalid purple segment: {}", toString());
        }
    }

    public String toString()
    {
        return String.format("loc(%s: %d - %d) ratio(%s) seg(%s)", Chromosome, mStart, mEnd, RatioSupport, Support);
    }
}
