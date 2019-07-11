package com.hartwig.hmftools.linx.types;

import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_PAIR;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.seIndex;

import java.util.List;

import com.google.common.collect.Lists;

public class SvChainState
{
    public final SvVarData SV;

    public final int MinPloidy;
    public final int Ploidy;
    public final int MaxPloidy;

    private int mCurrentCount;
    private int[] mBreakendCount;

    private final List<SvBreakend> mConnectionsStart;
    private final List<SvBreakend> mConnectionsEnd;

    final List<SvBreakend> mRepBreakendsStart;
    final List<SvBreakend> mRepBreakendsEnd;

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
            Ploidy = var.getImpliedPloidy();
            MinPloidy = (int) round(var.ploidyMin());
            MaxPloidy = (int) round(var.ploidyMax());
        }

        mCurrentCount = 0;
        mBreakendCount = new int[SE_PAIR];

        mConnectionsStart = Lists.newArrayList();
        mConnectionsEnd = Lists.newArrayList();
        mRepBreakendsEnd = Lists.newArrayList();
        mRepBreakendsStart = Lists.newArrayList();
    }

    public void add() { ++mCurrentCount; }
    public void add(boolean isStart) { ++mBreakendCount[seIndex(isStart)]; }

    public int curentCount() { return mCurrentCount; }
    public int breakendCount(int se) { return mBreakendCount[se]; }
    public int breakendCount(boolean isStart) { return mBreakendCount[seIndex(isStart)]; }

    public int minUnlinked() { return max(MinPloidy - mCurrentCount, 0); }
    public int maxUnlinked() { return max(MaxPloidy - mCurrentCount, 0); }
    public int unlinked() { return max(Ploidy - mCurrentCount, 0); }

    public int minUnlinked(boolean isStart) { return minUnlinked(seIndex(isStart)); }
    public int maxUnlinked(boolean isStart) { return maxUnlinked(seIndex(isStart)); }
    public int unlinked(boolean isStart) { return unlinked(seIndex(isStart)); }

    public int minUnlinked(int se) { return max(MinPloidy - mBreakendCount[se], 0); }
    public int maxUnlinked(int se) { return max(MaxPloidy - mBreakendCount[se],0); }
    public int unlinked(int se) { return max(Ploidy - mBreakendCount[se],0); }

    public final List<SvBreakend> getConnections(boolean isStart) { return isStart ? mConnectionsStart : mConnectionsEnd; }
    public int uniqueConnections(boolean isStart) { return isStart ? mConnectionsStart.size() : mConnectionsEnd.size(); }

    public void addConnection(final SvBreakend breakend, boolean isStart)
    {
        if(isStart && !mConnectionsStart.contains(breakend))
            mConnectionsStart.add(breakend);
        else if (!isStart && !mConnectionsEnd.contains(breakend))
            mConnectionsEnd.add(breakend);
    }

    public void addRepBreakend(final SvBreakend breakend)
    {
        if(breakend.usesStart() && !mRepBreakendsStart.contains(breakend))
            mRepBreakendsStart.add(breakend);
        else if (!breakend.usesStart() && !mRepBreakendsEnd.contains(breakend))
            mRepBreakendsEnd.add(breakend);
    }

    public void removeRepBreakend(final SvBreakend breakend)
    {
        if(breakend.usesStart())
            mRepBreakendsStart.remove(breakend);
        else if (!breakend.usesStart())
            mRepBreakendsEnd.remove(breakend);
    }

    public List<SvBreakend> getRepBreakends(boolean isStart) { return isStart ? mRepBreakendsStart : mRepBreakendsEnd; }

    public String toString()
    {
        return String.format("id(%d) ploidy(%d-%d-%d) counts(%d s=%d e=%d)",
                SV.dbId(), MinPloidy, Ploidy, MaxPloidy, mCurrentCount, mBreakendCount[SE_START], mBreakendCount[SE_END]);
    }

}
