package com.hartwig.hmftools.svanalysis.types;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.calcConsistency;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.makeChrArmStr;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SE_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isStart;

import com.google.common.collect.Lists;

import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// contains all the SVs for an arm within a cluster
public class SvArmGroup {

    private final String mId;
    private List<SvVarData> mSVs;
    private List<SvBreakend> mBreakends;
    final SvCluster mCluster;

    private final String mChromosome;
    private final String mArm;

    private boolean mRequiresRecalc;

    // width of SVs on the arm taking into account any excluded SVs
    private long mStartPos;
    private long mEndPos;
    private int mConsistency;

    private static final Logger LOGGER = LogManager.getLogger(SvArmGroup.class);

    public SvArmGroup(final SvCluster cluster, final String chr, final String arm)
    {
        mId = makeChrArmStr(chr, arm);
        mCluster = cluster;

        mChromosome = chr;
        mArm = arm;

        mSVs = Lists.newArrayList();
        mBreakends = Lists.newArrayList();
        mStartPos = -1;
        mEndPos = -1;
        mConsistency = 0;
        mRequiresRecalc = true;
    }

    public final String id() { return mId; }

    public final String chromosome() { return mChromosome; }
    public final String arm() { return mArm; }
    public long posStart() { return mStartPos; }
    public long posEnd() { return mEndPos; }

    public boolean hasEndsSet()
    {
        return mStartPos >= 0 && mEndPos >= 0;
    }

    public List<SvVarData> getSVs() { return mSVs; }
    public List<SvBreakend> getBreakends() { return mBreakends; }

    public void addVariant(final SvVarData var)
    {
        mSVs.add(var);

        for(int be = SE_START; be <= SE_END; ++be)
        {
            if(var.isNullBreakend() && be == SE_END)
                continue;

            SvBreakend breakend = var.getBreakend(isStart(be));

            if(!breakend.chromosome().equals(mChromosome) || !breakend.arm().equals(mArm))
                continue;

            int index = 0;
            while(index < mBreakends.size())
            {
                final SvBreakend otherBreakend = mBreakends.get(index);

                if(breakend.position() < otherBreakend.position())
                    break;

                ++index;
            }

            mBreakends.add(index, breakend);
        }

        mStartPos = mBreakends.get(0).position();
        mEndPos = mBreakends.get(mBreakends.size()-1).position();

        mRequiresRecalc = true;
    }

    public boolean matches(final SvArmGroup other)
    {
        return mChromosome.equals(other.chromosome()) && mArm.equals(other.arm());
    }

    public void setBoundaries(final List<SvVarData> bndIgnoreList)
    {
        if(!mRequiresRecalc)
            return;

        mConsistency = 0;

        for(final SvBreakend breakend : mBreakends)
        {
            mConsistency += calcConsistency(breakend.getSV(), breakend.usesStart());
        }

        //LOGGER.debug("cluster({}) arm({}) SVs({}) range({} -> {}) consistency({})",
        //        mCluster.id(), mId, mSVs.size(), mStartPos, mEndPos, mConsistency);

        mRequiresRecalc = false;
    }

    public boolean isConsistent() { return mConsistency == 0; }

    public boolean canLink(long maxBoundaryLength)
    {
        if(!isConsistent())
            return true;

        if(hasEndsSet() && posEnd() - posStart() >= maxBoundaryLength)
            return true;

        return false;
    }

}
