package com.hartwig.hmftools.svanalysis.types;

import static java.lang.Math.max;
import static java.lang.Math.min;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.svanalysis.analysis.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvClusterData;

import java.util.List;

public class SvArmGroup {

    private List<SvClusterData> mSVs;
    final SvCluster mCluster;

    private final String mChromosome;
    private final String mArm;
    private long mStartPos;
    private long mEndPos;

    public SvArmGroup(final SvCluster cluster, final String chr, final String arm)
    {
        mCluster = cluster;

        mChromosome = chr;
        mArm = arm;
        mStartPos = -1;
        mEndPos = -1;
        mSVs = Lists.newArrayList();
    }

    public final String posId() {
        return String.format("cl %d: %s_%s %d:%d", mCluster.getId(), mChromosome, mArm, mStartPos, mEndPos);
    }

    public final String chromosome() { return mChromosome; }
    public final String arm() { return mArm; }
    public long posStart() { return mStartPos; }
    public long posEnd() { return mEndPos; }

    public boolean hasEndsSet()
    {
        // for arms with only a BND or SGL
        return mStartPos >= 0 && mEndPos >= 0;
    }

    public List<SvClusterData> getSVs() { return mSVs; }
    public int getCount() { return mSVs.size(); }

    public void addVariant(final SvClusterData var)
    {
        mSVs.add(var);

        if(var.chromosome(true).equals(mChromosome))
        {
            mStartPos = mStartPos == -1 ? var.position(true) : min(mStartPos, var.position(true));
        }

        if(!var.isNullBreakend() && var.chromosome(false).equals(mChromosome))
        {
            mEndPos = max(mEndPos, var.position(false));
        }
    }

    public boolean matches(final SvArmGroup other)
    {
        return mChromosome.equals(other.chromosome()) && mArm.equals(other.arm());
    }
}
