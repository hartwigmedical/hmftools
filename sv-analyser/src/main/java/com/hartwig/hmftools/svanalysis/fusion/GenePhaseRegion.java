package com.hartwig.hmftools.svanalysis.fusion;

public class GenePhaseRegion
{
    public final String GeneId;
    public final int Phase;
    public final int RegionType;

    private long mStart;
    private long mEnd;

    public static final int REGION_TYPE_CODING = 0;
    public static final int REGION_TYPE_5PUTR = 1;
    public static final int REGION_TYPE_NON_CODING = 2;

    public GenePhaseRegion(final String geneId, final long start, final long end, final int phase, int regionType)
    {
        GeneId = geneId;
        Phase = phase;
        mStart = start;
        mEnd = end;
        RegionType = regionType;
    }

    public long start() { return mStart; }
    public void setStart(long start) { mStart = start; }

    public long end() { return mEnd; }
    public void setEnd(long end) { mEnd = end; }

    public long length() { return mEnd - mStart; }
}
