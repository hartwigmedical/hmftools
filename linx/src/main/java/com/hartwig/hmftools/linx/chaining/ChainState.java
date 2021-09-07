package com.hartwig.hmftools.linx.chaining;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.seIndex;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvVarData;

public class ChainState
{
    public final SvVarData SV;

    public final double MinJcn;
    public final double Jcn;
    public final double MaxJcn;

    private final double mExhaustedLevel;

    private final boolean mSingleJcn;

    private double[] mBreakendCount;

    // unique connections made to other SVs
    private final List<SvBreakend> mConnectionsStart;
    private final List<SvBreakend> mConnectionsEnd;

    public static final double EXHAUSTED_JCN_PERC = 0.1;
    public static final double EXHAUSTED_JCN_ABS = 0.2;

    public ChainState(final SvVarData var, Double singleJcn)
    {
        SV = var;

        if(singleJcn != null)
        {
            mSingleJcn = true;
            Jcn = singleJcn;
            MaxJcn = singleJcn;
            MinJcn = singleJcn;
        }
        else
        {
            mSingleJcn = false;
            Jcn = var.jcn();
            MinJcn = var.jcnMin();
            MaxJcn = max(var.jcnMax(), Jcn);
        }

        mExhaustedLevel = min(EXHAUSTED_JCN_PERC * Jcn, EXHAUSTED_JCN_ABS);

        mBreakendCount = new double[SE_PAIR];

        mConnectionsStart = Lists.newArrayList();
        mConnectionsEnd = Lists.newArrayList();
    }

    public void add(boolean isStart, double linkJcn) { mBreakendCount[seIndex(isStart)] += linkJcn; }
    public void set(boolean isStart, double linkJcn) { mBreakendCount[seIndex(isStart)] = linkJcn; }

    public double currentCount() { return min(mBreakendCount[SE_START], mBreakendCount[SE_END]); }
    public double breakendCount(boolean isStart) { return mBreakendCount[seIndex(isStart)]; }

    public boolean breakendExhausted(boolean isStart)
    {
        if(breakendCount(isStart) == 0)
            return false;

        return unlinked(isStart) <= mExhaustedLevel;
    }
    public boolean breakendExhaustedVsMax(boolean isStart) { return maxUnlinked(isStart) <= mExhaustedLevel; }

    public double unlinked() { return max(Jcn - currentCount(), 0); }
    public double unlinked(boolean isStart) { return unlinked(seIndex(isStart)); }
    public double unlinked(int se) { return max(Jcn - mBreakendCount[se],0); }

    public double maxUnlinked(boolean isStart) { return max(MaxJcn - mBreakendCount[seIndex(isStart)],0); }

    public final List<SvBreakend> getConnections(boolean isStart) { return isStart ? mConnectionsStart : mConnectionsEnd; }
    public int uniqueConnections(boolean isStart) { return isStart ? mConnectionsStart.size() : mConnectionsEnd.size(); }
    public boolean hasConnections() { return !mConnectionsStart.isEmpty() || !mConnectionsEnd.isEmpty(); }

    public void addConnection(final SvBreakend breakend, boolean isStart)
    {
        if(isStart && !mConnectionsStart.contains(breakend))
            mConnectionsStart.add(breakend);
        else if (!isStart && !mConnectionsEnd.contains(breakend))
            mConnectionsEnd.add(breakend);
    }

    public boolean overlaps(final ChainState other)
    {
        return !(MinJcn > other.MaxJcn || MaxJcn < other.MinJcn);
    }

    public boolean equals(final ChainState other)
    {
        return copyNumbersEqual(Jcn, other.Jcn);
    }

    public String toString()
    {
        if(mSingleJcn)
        {
            return String.format("id(%d) ploidy(%.1f) counts(s=%.1f e=%.1f)",
                    SV.id(), Jcn, mBreakendCount[SE_START], mBreakendCount[SE_END]);
        }
        else
        {
            if (Jcn > 10)
            {
                return String.format("id(%d) ploidy(%.0f-%.0f-%.0f) counts(s=%.0f e=%.0f)",
                        SV.id(), MinJcn, Jcn, MaxJcn, mBreakendCount[SE_START], mBreakendCount[SE_END]);
            }
            else if (Jcn < 0.5)
            {
                return String.format("id(%d) ploidy(%.2f-%.2f-%.2f) counts(s=%.2f e=%.2f)",
                        SV.id(), MinJcn, Jcn, MaxJcn, mBreakendCount[SE_START], mBreakendCount[SE_END]);
            }
            else
            {
                return String.format("id(%d) ploidy(%.1f-%.1f-%.1f) counts(s=%.1f e=%.1f)",
                        SV.id(), MinJcn, Jcn, MaxJcn, mBreakendCount[SE_START], mBreakendCount[SE_END]);
            }
        }
    }

}
