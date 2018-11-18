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

public class SvArmGroup {

    private final String mId;
    private List<SvVarData> mSVs;
    final SvCluster mCluster;

    private final String mChromosome;
    private final String mArm;

    // width of SVs on the arm taking into account any excluded SVs
    private boolean mRequiresRecalc;
    private long mStartPos;
    private long mEndPos;
    private int mConsistency;

    private SvBreakend mStartBreakend;
    private SvBreakend mEndBreakend;

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

        mStartBreakend = null;
        mEndBreakend = null;
    }

    public final String id() { return mId; }
    public final String posId() {
        return String.format("cl %d: %s_%s %d:%d", mCluster.getId(), mChromosome, mArm, mStartPos, mEndPos);
    }

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

    public void setBoundaries(final List<SvVarData> unlinkedBnds)
    {
        if(!mRequiresRecalc)
            return;

        mStartPos = -1;
        mEndPos = -1;
        mConsistency = 0;

        for(final SvVarData var : mSVs)
        {
            if(var.type() == BND && !unlinkedBnds.contains(var))
                continue;

            for(int be = SVI_START; be <= SVI_END; ++be)
            {
                if(be == SVI_END && var.isNullBreakend())
                    continue;

                boolean useStart = isStart(be);

                if(!var.chromosome(useStart).equals(mChromosome))
                    continue;

                mConsistency += calcConsistency(var, useStart);

                long position = var.position(useStart);

                mStartPos = mStartPos == -1 ? position : min(mStartPos, position);
                mEndPos = max(mEndPos, position);
            }
        }

        mRequiresRecalc = false;
    }

    public void setBreakend(SvBreakend breakend, boolean isStart)
    {
        if(isStart)
            mStartBreakend = breakend;
        else
            mEndBreakend = breakend;
    }

    public boolean noOpenBreakends()
    {
        return mStartBreakend == null && mEndBreakend == null;
    }

    public boolean isConsistent() { return mConsistency == 0; }

    /*
    public boolean isConsistent()
    {
        int consistency = calcConsistency(mSVs);

        if(consistency != 0)
            return false;

        // has both breakends set, face out towards telomere and centromere and have the same copy number
        // breakends can be from the same variant (should logically be a simple cluster-1)
        if(noOpenBreakends())
            return true;

        if(mStartBreakend == null || mEndBreakend == null)
            return false;

        if(!mStartBreakend.getSV().isSimpleType() && mStartBreakend.orientation() != 1)
            return false;

        if(!mEndBreakend.getSV().isSimpleType() && mEndBreakend.orientation() != -1)
            return false;

        double cnStart = mStartBreakend.getSV().copyNumber(mStartBreakend.usesStart());
        double cnEnd = mEndBreakend.getSV().copyNumber(mEndBreakend.usesStart());

        return copyNumbersEqual(cnStart, cnEnd);
    }
    */

}
