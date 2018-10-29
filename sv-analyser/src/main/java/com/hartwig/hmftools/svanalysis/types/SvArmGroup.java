package com.hartwig.hmftools.svanalysis.types;

import static java.lang.Math.max;
import static java.lang.Math.min;

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
    private long mStartPos;
    private long mEndPos;

    private SvBreakend mStartBreakend;
    private SvBreakend mEndBreakend;

    public SvArmGroup(final SvCluster cluster, final String chr, final String arm)
    {
        mId = makeChrArmStr(chr, arm);
        mCluster = cluster;

        mChromosome = chr;
        mArm = arm;
        mStartPos = -1;
        mEndPos = -1;
        mSVs = Lists.newArrayList();

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
        // for arms with only a BND or SGL
        return mStartPos >= 0 && mEndPos >= 0;
    }

    public List<SvVarData> getSVs() { return mSVs; }
    public int getCount() { return mSVs.size(); }

    public void addVariant(final SvVarData var)
    {
        mSVs.add(var);

        for(int be = SVI_START; be <= SVI_END; ++be)
        {
            if(be == SVI_END && var.isNullBreakend())
                continue;

            boolean useStart = isStart(be);

            if(!var.chromosome(useStart).equals(mChromosome))
                continue;

            long position = var.position(useStart);

            mStartPos = mStartPos == -1 ? position : min(mStartPos, position);
            mEndPos = max(mEndPos, position);
        }
    }

    public boolean matches(final SvArmGroup other)
    {
        return mChromosome.equals(other.chromosome()) && mArm.equals(other.arm());
    }

    public SvBreakend getBreakend(boolean useStart) { return useStart ? mStartBreakend : mEndBreakend; }

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

    public boolean isConsistent()
    {
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

}
