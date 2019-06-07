package com.hartwig.hmftools.svanalysis.fusion;

public class GenePhaseRegion
{
    public final String GeneId;
    public final int Phase;

    private long mStart;
    private long mEnd;

    public GenePhaseRegion(final String geneId, final long start, final long end, final int phase)
    {
        GeneId = geneId;
        Phase = phase;
        mStart = start;
        mEnd = end;
    }

    public long start() { return mStart; }
    public void setStart(long start) { mStart = start; }

    public long end() { return mEnd; }
    public void setEnd(long end) { mEnd = end; }
}
