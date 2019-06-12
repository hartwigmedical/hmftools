package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.linx.fusion.GenePhaseType.PHASE_0;
import static com.hartwig.hmftools.linx.fusion.GenePhaseType.PHASE_1;
import static com.hartwig.hmftools.linx.fusion.GenePhaseType.PHASE_2;
import static com.hartwig.hmftools.linx.fusion.GenePhaseType.PHASE_5P_UTR;
import static com.hartwig.hmftools.linx.fusion.GenePhaseType.PHASE_MAX;
import static com.hartwig.hmftools.linx.fusion.GenePhaseType.PHASE_NON_CODING;
import static com.hartwig.hmftools.linx.fusion.GenePhaseType.typeAsInt;

public class GenePhaseRegion
{
    public final String GeneId;
    public GenePhaseType Phase;
    public final int RegionType;

    private long mStart;
    private long mEnd;

    private boolean[] mPhaseArray;
    private boolean[] mPreGenePhaseStatus;
    private int mCombinedPhase;

    public static final int REGION_TYPE_CODING = 0;
    public static final int REGION_TYPE_5PUTR = 1;
    public static final int REGION_TYPE_NON_CODING = 2;
    public static final int REGION_TYPE_MIXED = 3;

    public GenePhaseRegion(final String geneId, long start, long end, GenePhaseType phase)
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

        mPhaseArray = new boolean[GenePhaseType.values().length];
        mPreGenePhaseStatus = new boolean[GenePhaseType.values().length];

        mPhaseArray[typeAsInt(phase)] = true;
        calcCombinedPhase();
    }

    public GenePhaseRegion(final String geneId, final long start, final long end, final boolean[] phaseArray, final boolean[] preGeneArray)
    {
        GeneId = geneId;
        Phase = PHASE_5P_UTR; // will not be used
        mStart = start;
        mEnd = end;

        RegionType = REGION_TYPE_MIXED;

        mPhaseArray = new boolean[PHASE_MAX];
        mPreGenePhaseStatus = new boolean[GenePhaseType.values().length];
        addPhases(phaseArray, preGeneArray);
        calcCombinedPhase();
    }

    public void setPreGene(boolean toggle, GenePhaseType phase)
    {
        mPreGenePhaseStatus[typeAsInt(phase)] = true;
    }

    public static GenePhaseType mapExonPhase(int exonPhase)
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
    public final boolean[] getPreGenePhaseStatus() { return mPreGenePhaseStatus; }

    public void addPhases(final boolean[] phases, final boolean[] preGeneArray)
    {
        if(phases.length != mPhaseArray.length)
            return;

        for(int i = 0; i < PHASE_MAX; ++i)
        {
            mPhaseArray[i] |= phases[i];
            mPreGenePhaseStatus[i] |= preGeneArray[i];
        }

        calcCombinedPhase();
    }

    public boolean hasPhase(GenePhaseType phase)
    {
        return mPhaseArray[typeAsInt(phase)];
    }

    public boolean hasPhaseOnly(GenePhaseType phase)
    {
        int phaseInt = typeAsInt(phase);
        for(int i = 0; i < PHASE_MAX; ++i)
        {
            if(mPhaseArray[i] && phaseInt != i)
                return false;
            if(!mPhaseArray[i] && phaseInt == i)
                return false;
        }

        return true;
    }

    public static boolean hasAnyPhaseMatch(final GenePhaseRegion regionUp, final GenePhaseRegion regionDown, boolean allowPreGeneDown)
    {
        for(int i = 0; i < PHASE_MAX; ++i)
        {
            if(regionUp.getPhaseArray()[i] && regionDown.getPhaseArray()[i])
            {
                if((regionUp.getPreGenePhaseStatus()[i]))
                    continue;

                if((!allowPreGeneDown && regionDown.getPreGenePhaseStatus()[i]))
                    continue;

                return true;
            }
        }

        return false;
    }

    public static  boolean regionsPhaseMatched(final GenePhaseRegion regionUp, final GenePhaseRegion regionDown)
    {
        if(regionUp.hasPhase(PHASE_NON_CODING) && regionDown.hasPhase(PHASE_5P_UTR))
            return true;

        if(hasAnyPhaseMatch(regionUp, regionDown, true))
            return true;

        return false;
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

    public boolean isAnyPreGene()
    {
        for(int i = 0; i < PHASE_MAX; ++i)
        {
            if(mPreGenePhaseStatus[i])
                return true;
        }

        return false;
    }

    public boolean hasNonPreGene()
    {
        for(int i = 0; i < PHASE_MAX; ++i)
        {
            if(mPhaseArray[i] && !mPreGenePhaseStatus[i])
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

    public static int simpleToCombinedPhase(GenePhaseType phase)
    {
        return (int)Math.pow(10, typeAsInt(phase));

    }

}
