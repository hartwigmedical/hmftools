package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.linx.fusion.GenePhaseType.PHASE_0;
import static com.hartwig.hmftools.linx.fusion.GenePhaseType.PHASE_1;
import static com.hartwig.hmftools.linx.fusion.GenePhaseType.PHASE_2;
import static com.hartwig.hmftools.linx.fusion.GenePhaseType.PHASE_5P_UTR;
import static com.hartwig.hmftools.linx.fusion.GenePhaseType.PHASE_MAX;
import static com.hartwig.hmftools.linx.fusion.GenePhaseType.PHASE_NON_CODING;
import static com.hartwig.hmftools.linx.fusion.GenePhaseType.intAsType;
import static com.hartwig.hmftools.linx.fusion.GenePhaseType.typeAsInt;

import java.util.StringJoiner;

public class GenePhaseRegion
{
    public final String GeneId;
    public GenePhaseType Phase;

    private long mStart;
    private long mEnd;

    private boolean[] mPhaseArray;
    private boolean[] mPreGenePhaseStatus;
    private int mCombinedPhase;

    public GenePhaseRegion(final String geneId, long start, long end, GenePhaseType phase)
    {
        GeneId = geneId;
        Phase = phase;
        mStart = start;
        mEnd = end;

        mPhaseArray = new boolean[GenePhaseType.values().length];
        mPreGenePhaseStatus = new boolean[GenePhaseType.values().length];

        mPhaseArray[typeAsInt(phase)] = true;
        calcCombinedPhase();
    }

    public GenePhaseRegion(final String geneId, long start, long end, final boolean[] phaseArray, final boolean[] preGeneArray)
    {
        GeneId = geneId;
        Phase = PHASE_5P_UTR; // will not be used
        mStart = start;
        mEnd = end;

        mPhaseArray = new boolean[PHASE_MAX];
        mPreGenePhaseStatus = new boolean[GenePhaseType.values().length];
        addPhases(phaseArray, preGeneArray);
        calcCombinedPhase();
    }

    public void setPreGene(boolean toggle, GenePhaseType phase)
    {
        mPreGenePhaseStatus[typeAsInt(phase)] = toggle;
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

    public static final String PD_DELIMITER = ":";

    private static String boolToStr(boolean value) { return value ? "1" : "0"; }
    private static boolean strToBool(final String value) { return value.equals("1"); }

    public String toCsv(boolean useArray)
    {
        StringJoiner output = new StringJoiner(PD_DELIMITER);
        output.add(String.valueOf(mStart));
        output.add(String.valueOf(mEnd));

        if(useArray)
        {
            for (int i = 0; i < PHASE_MAX; ++i)
            {
                output.add(boolToStr(mPhaseArray[i]));
            }

            for (int i = 0; i < PHASE_MAX; ++i)
            {
                output.add(boolToStr(mPreGenePhaseStatus[i]));
            }
        }
        else
        {
            int phase = typeAsInt(Phase);
            output.add(String.valueOf(phase));
            output.add(boolToStr(mPreGenePhaseStatus[phase]));
        }

        return output.toString();
    }

    public static GenePhaseRegion fromCsv(final String geneId, final String inputStr, boolean useArray)
    {
        String[] items = inputStr.split(PD_DELIMITER);

        int startEndItems = 2;

        if(useArray)
        {
            if (items.length != startEndItems + PHASE_MAX * 2)
                return null;
        }
        else
        {
            if(items.length != 4)
                return null;
        }

        long start = Long.parseLong(items[0]);
        long end = Long.parseLong(items[1]);

        if(useArray)
        {
            boolean[] phases = new boolean[PHASE_MAX];
            boolean[] status = new boolean[PHASE_MAX];

            for (int i = 0; i < PHASE_MAX; ++i)
            {
                phases[i] = strToBool(items[i + startEndItems]);
                status[i] = strToBool(items[i + startEndItems + PHASE_MAX]);
            }

            return new GenePhaseRegion(geneId, start, end, phases, status);
        }
        else
        {
            GenePhaseType phase = intAsType(Integer.parseInt(items[2]));
            boolean preGene = strToBool(items[3]);

            GenePhaseRegion region = new GenePhaseRegion(geneId, start, end, phase);
            region.setPreGene(preGene, phase);
            return region;
        }
    }

}
