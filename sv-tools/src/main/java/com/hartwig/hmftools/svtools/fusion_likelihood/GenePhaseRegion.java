package com.hartwig.hmftools.svtools.fusion_likelihood;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.svtools.fusion_likelihood.GenePhaseType.PHASE_0;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GenePhaseType.PHASE_1;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GenePhaseType.PHASE_2;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GenePhaseType.PHASE_5P_UTR;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GenePhaseType.PHASE_MAX;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GenePhaseType.PHASE_NON_CODING;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GenePhaseType.typeAsInt;

import java.util.List;

public class GenePhaseRegion
{
    public final String GeneId;
    public GenePhaseType Phase;

    private int mStart;
    private int mEnd;

    private boolean[] mPhaseArray;
    private boolean[] mPreGenePhaseStatus;
    private int mCombinedPhase;
    private int mCombinedPreGeneStatus;

    private boolean mHasOverlaps;
    private int mTransId;
    private boolean mProteinCoding;

    public GenePhaseRegion(final String geneId, int start, int end, GenePhaseType phase)
    {
        GeneId = geneId;
        Phase = phase;
        mStart = start;
        mEnd = end;

        mHasOverlaps = false;
        mTransId = 0;
        mProteinCoding = false;

        mPhaseArray = new boolean[GenePhaseType.values().length];
        mPreGenePhaseStatus = new boolean[GenePhaseType.values().length];

        mPhaseArray[typeAsInt(phase)] = true;
        calcCombinedPhase();
    }

    public GenePhaseRegion(final String geneId, int start, int end, final boolean[] phaseArray, final boolean[] preGeneArray)
    {
        GeneId = geneId;
        Phase = PHASE_5P_UTR; // will not be used
        mStart = start;
        mEnd = end;
        mTransId = 0;

        mPhaseArray = new boolean[PHASE_MAX];
        mPreGenePhaseStatus = new boolean[GenePhaseType.values().length];
        addPhases(phaseArray, preGeneArray);
        calcCombinedPhase();
    }

    public static GenePhaseRegion from(final GenePhaseRegion other, int start, int end)
    {
        // copy with new region
        GenePhaseRegion newRegion = new GenePhaseRegion(
                other.GeneId, start, end, other.getPhaseArray(), other.getPreGenePhaseStatus());

        newRegion.setProteinCoding(other.proteinCoding());

        return newRegion;
    }

    public static GenePhaseRegion from(final GenePhaseRegion other)
    {
        // exact copy
        return GenePhaseRegion.from(other, other.start(), other.end());
    }

    public void setPreGene(boolean toggle, GenePhaseType phase)
    {
        mPreGenePhaseStatus[typeAsInt(phase)] = toggle;
        calcCombinedPhase();
    }

    public boolean hasPreGeneStatus() { return mCombinedPreGeneStatus > 0; }

    public void setProteinCoding(boolean toggle) { mProteinCoding = toggle; }
    public boolean proteinCoding() { return mProteinCoding; }

    public void setHasOverlaps(boolean toggle) { mHasOverlaps = toggle; }
    public boolean hasOverlaps() { return mHasOverlaps; }

    public void setTransId(int transId) { mTransId = transId; }
    public int transId() { return mTransId; }

    public static GenePhaseType mapExonPhase(int exonPhase)
    {
        if (exonPhase == 0)
            return PHASE_0;
        else if (exonPhase == 1)
            return PHASE_1;
        else if (exonPhase == 2)
            return PHASE_2;
        else
            return PHASE_5P_UTR;
    }

    public int start() { return mStart; }
    public void setStart(int start) { mStart = start; }

    public int end() { return mEnd; }
    public void setEnd(int end) { mEnd = end; }

    public int length() { return mEnd - mStart; }

    public final boolean[] getPhaseArray() { return mPhaseArray; }
    public final boolean[] getPreGenePhaseStatus() { return mPreGenePhaseStatus; }

    public void addPhases(final boolean[] phases, final boolean[] preGeneArray)
    {
        if (phases.length != mPhaseArray.length)
            return;

        for (int i = 0; i < PHASE_MAX; ++i)
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
        // look for any phase match but exclude non-coding since this is handled in the direction-specific method below
        for(int i = 0; i < PHASE_MAX; ++i)
        {
            if(i == typeAsInt(PHASE_NON_CODING))
                continue;

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

    public static boolean regionsPhaseMatched(final GenePhaseRegion regionUp, final GenePhaseRegion regionDown)
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

    public boolean hasPhasedType()
    {
        return hasPhase(PHASE_0) || hasPhase(PHASE_1) || hasPhase(PHASE_2);
    }

    public int getCombinedPhase() { return mCombinedPhase; }

    private void calcCombinedPhase()
    {
        mCombinedPhase = calcCombinedPhase(mPhaseArray);
        mCombinedPreGeneStatus = calcCombinedPhase(mPreGenePhaseStatus);
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

    public int getCombinedPreGeneStatus() { return mCombinedPreGeneStatus; }

    public static int simpleToCombinedPhase(GenePhaseType phase)
    {
        return (int)Math.pow(10, typeAsInt(phase));
    }

    public void populateLengthCounts(int[] counts, boolean allowPreGene)
    {
        if(counts.length != PHASE_MAX)
            return;

        for(int i = 0; i < PHASE_MAX; ++i)
        {
            if(mPhaseArray[i] && allowPreGene == mPreGenePhaseStatus[i])
                counts[i] += length();
        }
    }

    public static boolean haveOverlap(final GenePhaseRegion region1, final GenePhaseRegion region2, int overlapBases)
    {
        // a positive value for overlap bases requires a gap between the regions, a negative value allows some overlap without calling it such
        if (region1.end() < region2.start() - overlapBases || region1.start() > region2.end() + overlapBases)
            return false;

        return true;
    }

    public static boolean hasNoOverlappingRegions(final List<GenePhaseRegion> regions)
    {
        for (int i = 0; i < regions.size(); ++i)
        {
            GenePhaseRegion region1 = regions.get(i);

            for (int j = i + 1; j < regions.size(); ++j)
            {
                GenePhaseRegion region2 = regions.get(j);

                if(haveOverlap(region1, region2, 0))
                    return false;
            }
        }

        return true;
    }

    public final String toString()
    {
        return String.format("%s: range(%d - %d) len(%d) phases(%d) preGene(%d)",
                GeneId, mStart, mEnd, length(), mCombinedPhase, mCombinedPreGeneStatus);
    }
}
