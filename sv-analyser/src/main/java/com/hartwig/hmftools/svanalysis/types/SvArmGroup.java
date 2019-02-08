package com.hartwig.hmftools.svanalysis.types;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.calcConsistency;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.makeChrArmStr;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isStart;

import com.google.common.collect.Lists;

import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// contains all the SVs for an arm within a cluster
public class SvArmGroup {

    private final String mId;
    private List<SvVarData> mSVs;
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

    public void addVariant(final SvVarData var)
    {
        mSVs.add(var);
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

        mStartPos = -1;
        mEndPos = -1;
        mConsistency = 0;

        for(final SvVarData var : mSVs)
        {
            if(var.type() == BND && bndIgnoreList.contains(var)) // skip any BNDs only in this arm as a short TI
                continue;

            for(int be = SVI_START; be <= SVI_END; ++be)
            {
                if(be == SVI_END && var.isNullBreakend())
                    continue;

                boolean useStart = isStart(be);

                if(!var.chromosome(useStart).equals(mChromosome) || !var.arm(useStart).equals(mArm))
                    continue;

                mConsistency += calcConsistency(var, useStart);

                long position = var.position(useStart);

                mStartPos = mStartPos == -1 ? position : min(mStartPos, position);
                mEndPos = max(mEndPos, position);
            }
        }

        LOGGER.debug("cluster({}) arm({}) SVs({}) range({} -> {}) consistency({})",
                mCluster.id(), mId, mSVs.size(), mStartPos, mEndPos, mConsistency);

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
