package com.hartwig.hmftools.linx.chaining;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_PAIR;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.seIndex;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvVarData;

public class SvChainState
{
    public final SvVarData SV;

    public final double MinPloidy;
    public final double Ploidy;
    public final double MaxPloidy;

    private final double mExhaustedLevel;

    private double[] mBreakendCount;

    // unique connections made to other SVs
    private final List<SvBreakend> mConnectionsStart;
    private final List<SvBreakend> mConnectionsEnd;

    private static final double EXHAUSTED_PLOIDY_PERC = 0.1;
    private static final double EXHAUSTED_PLOIDY_ABS = 0.2;

    public SvChainState(final SvVarData var, boolean singlePloidy)
    {
        SV = var;

        if(singlePloidy)
        {
            Ploidy = 1;
            MaxPloidy = 1;
            MinPloidy = 1;
        }
        else
        {
            Ploidy = var.ploidy();
            MinPloidy = var.ploidyMin();
            MaxPloidy = max(var.ploidyMax(), Ploidy);
        }

        mExhaustedLevel = min(EXHAUSTED_PLOIDY_PERC * Ploidy, EXHAUSTED_PLOIDY_ABS);

        mBreakendCount = new double[SE_PAIR];

        mConnectionsStart = Lists.newArrayList();
        mConnectionsEnd = Lists.newArrayList();
    }

    public void add(boolean isStart, double linkPloidy) { mBreakendCount[seIndex(isStart)] += linkPloidy; }

    public double curentCount() { return min(mBreakendCount[SE_START], mBreakendCount[SE_END]); }
    public double breakendCount(boolean isStart) { return mBreakendCount[seIndex(isStart)]; }

    public boolean breakendExhausted(boolean isStart)
    {
        if(breakendCount(isStart) == 0)
            return false;

        return unlinked(isStart) <= mExhaustedLevel;
    }
    public boolean breakendExhaustedVsMax(boolean isStart) { return maxUnlinked(isStart) <= mExhaustedLevel; }

    public double unlinked() { return max(Ploidy - curentCount(), 0); }
    public double unlinked(boolean isStart) { return unlinked(seIndex(isStart)); }
    public double unlinked(int se) { return max(Ploidy - mBreakendCount[se],0); }

    public double minUnlinked() { return max(MinPloidy - curentCount(), 0); }
    public double maxUnlinked() { return max(MaxPloidy - curentCount(), 0); }
    public double minUnlinked(int se) { return max(MinPloidy - mBreakendCount[se], 0); }
    public double maxUnlinked(boolean isStart) { return max(MaxPloidy - mBreakendCount[seIndex(isStart)],0); }

    public final List<SvBreakend> getConnections(boolean isStart) { return isStart ? mConnectionsStart : mConnectionsEnd; }
    public int uniqueConnections(boolean isStart) { return isStart ? mConnectionsStart.size() : mConnectionsEnd.size(); }

    public void addConnection(final SvBreakend breakend, boolean isStart)
    {
        if(isStart && !mConnectionsStart.contains(breakend))
            mConnectionsStart.add(breakend);
        else if (!isStart && !mConnectionsEnd.contains(breakend))
            mConnectionsEnd.add(breakend);
    }

    public boolean overlaps(final SvChainState other)
    {
        return !(MinPloidy > other.MaxPloidy || MaxPloidy < other.MinPloidy);
    }

    public boolean equals(final SvChainState other)
    {
        return copyNumbersEqual(Ploidy, other.Ploidy);
    }

    public String toString()
    {
        if(Ploidy > 10)
        {
            return String.format("id(%d) ploidy(%.0f-%.0f-%.0f) counts(s=%.0f e=%.0f)",
                    SV.dbId(), MinPloidy, Ploidy, MaxPloidy, mBreakendCount[SE_START], mBreakendCount[SE_END]);
        }
        else if(Ploidy < 0.5)
        {
            return String.format("id(%d) ploidy(%.2f-%.2f-%.2f) counts(s=%.2f e=%.2f)",
                    SV.dbId(), MinPloidy, Ploidy, MaxPloidy, mBreakendCount[SE_START], mBreakendCount[SE_END]);
        }
        else
        {
            return String.format("id(%d) ploidy(%.1f-%.1f-%.1f) counts(s=%.1f e=%.1f)",
                    SV.dbId(), MinPloidy, Ploidy, MaxPloidy, mBreakendCount[SE_START], mBreakendCount[SE_END]);
        }
    }

}
