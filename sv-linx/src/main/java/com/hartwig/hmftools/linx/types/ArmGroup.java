package com.hartwig.hmftools.linx.types;

import static com.hartwig.hmftools.linx.analysis.SvUtilities.makeChrArmStr;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.segment.ChromosomeArm;

import java.util.List;

// contains all the SVs for an arm within a cluster
public class ArmGroup
{
    private final String mId;
    private List<SvVarData> mSVs;
    private List<SvBreakend> mBreakends;

    private final String mChromosome;
    private final ChromosomeArm mArm;

    // width of SVs on the arm taking into account any excluded SVs
    private int mStartPos;
    private int mEndPos;

    public ArmGroup(final String chr, final ChromosomeArm arm)
    {
        mId = makeChrArmStr(chr, arm);

        mChromosome = chr;
        mArm = arm;

        mSVs = Lists.newArrayList();
        mBreakends = Lists.newArrayList();
        mStartPos = -1;
        mEndPos = -1;
    }

    public final String id() { return mId; }

    public final String chromosome() { return mChromosome; }
    public final ChromosomeArm arm() { return mArm; }
    public int posStart() { return mStartPos; }
    public int posEnd() { return mEndPos; }

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
            if(var.isSglBreakend() && be == SE_END)
                continue;

            SvBreakend breakend = var.getBreakend(be);

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
    }

    public boolean matches(final ArmGroup other)
    {
        return mChromosome.equals(other.chromosome()) && mArm.equals(other.arm());
    }

}
