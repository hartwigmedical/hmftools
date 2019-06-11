package com.hartwig.hmftools.svanalysis.fusion;

public class GenePhaseRegion
{
    public final String GeneId;
    public int Phase;
    public final int RegionType;

    private long mStart;
    private long mEnd;

    private boolean[] mPhaseArray;
    private int mCombinedPhase;

    public static final int REGION_TYPE_CODING = 0;
    public static final int REGION_TYPE_5PUTR = 1;
    public static final int REGION_TYPE_NON_CODING = 2;
    public static final int REGION_TYPE_MIXED = 3;

    public static final int PHASE_NON_CODING = 0;
    public static final int PHASE_5P_UTR = 1;
    public static final int PHASE_0 = 2;
    public static final int PHASE_1 = 3;
    public static final int PHASE_2 = 4;
    public static final int PHASE_MAX = 5;

    public GenePhaseRegion(final String geneId, long start, long end, int phase)
    {
        GeneId = geneId;
        Phase = phase;
        mStart = start;
        mEnd = end;

        if(phase == PHASE_NON_CODING)
            RegionType = REGION_TYPE_NON_CODING;
        else if(phase == PHASE_5P_UTR)
            RegionType = REGION_TYPE_5PUTR;
        else
            RegionType = REGION_TYPE_CODING;

        mPhaseArray = new boolean[PHASE_MAX];

        if(!validPhase(phase))
            return;

        mPhaseArray[phase] = true;
        calcCombinedPhase();
    }

    public GenePhaseRegion(final String geneId, final long start, final long end, final boolean[] phaseArray)
    {
        GeneId = geneId;
        Phase = 0; // will not be used
        mStart = start;
        mEnd = end;

        RegionType = REGION_TYPE_MIXED;

        mPhaseArray = new boolean[PHASE_MAX];
        addPhases(phaseArray);
        calcCombinedPhase();
    }

    public static int mapExonPhase(int exonPhase)
    {
        if(exonPhase == 0)
            return PHASE_0;
        else if(exonPhase == 1)
            return PHASE_1;
        else if(exonPhase == 2)
            return PHASE_2;
        else
            return PHASE_5P_UTR;
    }

    public long start() { return mStart; }
    public void setStart(long start) { mStart = start; }

    public long end() { return mEnd; }
    public void setEnd(long end) { mEnd = end; }

    public long length() { return mEnd - mStart; }

    public final boolean[] getPhaseArray() { return mPhaseArray; }

    public void addPhases(final boolean[] phases)
    {
        if(phases.length != mPhaseArray.length)
            return;

        for(int i = 0; i < PHASE_MAX; ++i)
        {
            mPhaseArray[i] |= phases[i];
        }

        calcCombinedPhase();
    }

    public static boolean validPhase(int phase) { return phase >= 0 && phase < PHASE_MAX; }

    public boolean hasPhase(int phase)
    {
        if(!validPhase(phase))
            return false;

        return mPhaseArray[phase];
    }

    public boolean hasPhaseOnly(int phase)
    {
        if(!validPhase(phase))
            return false;

        for(int i = 0; i < PHASE_MAX; ++i)
        {
            if(mPhaseArray[i] && phase != i)
                return false;
            if(!mPhaseArray[i] && phase == i)
                return false;
        }

        return true;
    }

    public boolean hasAnyPhaseMatch(final boolean[] phaseArray)
    {
        for(int i = 0; i < PHASE_MAX; ++i)
        {
            if(mPhaseArray[i] && phaseArray[i])
                return true;
        }

        return false;
    }


    public int getCombinedPhase() { return mCombinedPhase; }

    private void calcCombinedPhase()
    {
        mCombinedPhase = calcCombinedPhase(mPhaseArray);
    }

    public static int calcCombinedPhase(final boolean[] phases)
    {
        if(phases.length != PHASE_MAX)
            return -1;

        int combinedPhase = 0;
        for(int i = 0; i < PHASE_MAX; ++i)
        {
            if(phases[i])
                combinedPhase += Math.pow(10, i);
        }

        return combinedPhase;
    }

    public static int simpleToCombinedPhase(int phase)
    {
        if(!validPhase(phase))
            return -1;

        return (int)Math.pow(10, phase);

    }

}
