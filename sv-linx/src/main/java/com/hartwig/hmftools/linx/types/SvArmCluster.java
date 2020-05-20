package com.hartwig.hmftools.linx.types;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.makeChrArmStr;
import static com.hartwig.hmftools.linx.chaining.LinkFinder.getMinTemplatedInsertionLength;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// clusters of proximate SVs on an arm
public class SvArmCluster
{
    private int mId;
    private List<SvBreakend> mBreakends;
    final SvCluster mCluster;

    private final String mChromosome;
    private final ChromosomeArm mArm;
    private int mStartPos;
    private int mEndPos;

    private int mType;
    private int mTICount;

    private static final Logger LOGGER = LogManager.getLogger(SvArmCluster.class);

    public SvArmCluster(int id, final SvCluster cluster, final String chr, final ChromosomeArm arm)
    {
        mId = id;
        mCluster = cluster;
        mArm = arm;
        mChromosome = chr;
        mStartPos = 0;
        mEndPos = 0;
        mBreakends = Lists.newArrayList();
        mType = -1;
        mTICount = 0;
    }

    public int id() { return mId; }

    public String getChrArm() { return makeChrArmStr(mChromosome, mArm); }

    public final String toString()
    {
        return String.format("%d: %s %s_%s %d:%d",
                mId, typeToString(mType), mChromosome, mArm, mStartPos, mEndPos);
    }

    public final String chromosome() { return mChromosome; }
    public final ChromosomeArm arm() { return mArm; }
    public int posStart() { return mStartPos; }
    public int posEnd() { return mEndPos; }

    public List<SvBreakend> getBreakends() { return mBreakends; }

    public void addBreakend(final SvBreakend breakend)
    {
        int index = 0;
        while(index < mBreakends.size())
        {
            final SvBreakend be = mBreakends.get(index);

            if(breakend.position() < be.position())
                break;

            ++index;
        }

        mBreakends.add(index, breakend);

        mStartPos = mBreakends.get(0).position();
        mEndPos = mBreakends.get(mBreakends.size()-1).position();
    }

    public static final int ARM_CL_ISOLATED_BE = 0;
    public static final int ARM_CL_TI_ONLY = 1;
    public static final int ARM_CL_DSB = 2; // includes simple DELs within the arm-cluster window
    public static final int ARM_CL_FOLDBACK = 3;
    public static final int ARM_CL_FOLDBACK_DSB = 4;
    public static final int ARM_CL_COMPLEX_FOLDBACK = 5;
    public static final int ARM_CL_COMPLEX_LINE = 6;
    public static final int ARM_CL_COMPLEX_OTHER = 7;
    public static final int ARM_CL_SIMPLE_DUP = 8;
    public static final int ARM_CL_SAME_ORIENT = 9;
    public static final int ARM_CL_MAX = ARM_CL_SAME_ORIENT;

    public static final String typeToString(int type)
    {
        switch(type)
        {
            case ARM_CL_ISOLATED_BE: return "ISOLATED_BE";
            case ARM_CL_TI_ONLY: return "TI_ONLY";
            case ARM_CL_DSB : return "DSB";
            case ARM_CL_FOLDBACK : return "FOLDBACK";
            case ARM_CL_FOLDBACK_DSB : return "FOLDBACK_DSB";
            case ARM_CL_COMPLEX_FOLDBACK : return "COMPLEX_FOLDBACK";
            case ARM_CL_COMPLEX_LINE: return "COMPLEX_LINE";
            case ARM_CL_COMPLEX_OTHER : return "COMPLEX_OTHER";
            case ARM_CL_SIMPLE_DUP: return "SIMPLE_DUP";
            case ARM_CL_SAME_ORIENT: return "SAME_ORIENT";
        }

        return "UNKNOWN";
    }

    public int getType() { return mType; }
    public String getTypeStr() { return typeToString(mType); }
    public int getTICount() { return mTICount; }

    public void setFeatures()
    {
        mTICount = 0;

        if(mBreakends.size() == 1)
        {
            mType = ARM_CL_ISOLATED_BE;
            return;
        }

        if(mBreakends.size() == 2)
        {
            final SvBreakend be1 = mBreakends.get(0);
            final SvVarData var1 = be1.getSV();
            final SvBreakend be2 = mBreakends.get(1);
            final SvVarData var2 = be2.getSV();

            if(var1 == var2)
            {
                if(var1.type() == DEL)
                {
                    mType = ARM_CL_DSB;
                    return;
                }
                else if(var1.type() == DUP)
                {
                    mType = ARM_CL_SIMPLE_DUP;
                    return;
                }
            }

            if(var1.getFoldbackId(be1.usesStart()) == var2.id())
            {
                mType = ARM_CL_FOLDBACK;
                return;
            }

            final SvLinkedPair tiPair = var1.getLinkedPair(be1.usesStart());

            if(tiPair != null && tiPair.hasBreakend(be2))
            {
                mType = ARM_CL_TI_ONLY;
                mTICount = 1;
                return;
            }

            final SvLinkedPair dbPair = be1.getSV().getDBLink(be1.usesStart());

            if(dbPair != null && dbPair == be2.getSV().getDBLink(be2.usesStart()))
            {
                mType = ARM_CL_DSB;
                return;
            }
        }

        // otherwise count the number of foldbacks, LINE, TIs and DSBs to determine the type
        int dsbCount = 0;
        int foldbackCount = 0;
        int suspectLine = 0;
        int tiCount = 0;

        List<SvBreakend> tiBreakends = Lists.newArrayList();

        for (final SvBreakend breakend : mBreakends)
        {
            if(breakend.getSV().isLineElement(breakend.usesStart()))
            {
                ++suspectLine;
            }

            if(breakend.getSV().isFoldback())
            {
                ++foldbackCount;
            }

            if(breakend.orientation() == -1)
            {
                final SvLinkedPair tiPair = breakend.getSV().getLinkedPair(breakend.usesStart());

                if (tiPair != null && mBreakends.contains(tiPair.getOtherBreakend(breakend)))
                {
                    tiBreakends.add(breakend);
                    tiBreakends.add(tiPair.getOtherBreakend(breakend));
                    ++tiCount;
                }
            }
        }

        if((foldbackCount % 2) == 2)
        {
            LOGGER.warn("cluster({}) armGroup({}) has invalid foldback count({})", mCluster.id(), toString(), foldbackCount);
        }
        else
        {
            foldbackCount /= 2; // to avoid double counting
        }

        mTICount = tiCount;

        if(tiBreakends.size() == mBreakends.size())
        {
            mType = ARM_CL_TI_ONLY;
            return;
        }

        if(suspectLine > 0)
        {
            mType = ARM_CL_COMPLEX_LINE;
            return;
        }
        else if(foldbackCount > 1)
        {
            mType = ARM_CL_COMPLEX_FOLDBACK;
            return;
        }

        int nonTiBreakendCount = mBreakends.size() - tiCount * 2;

        /*
        if(foldbackCount == 0 && nonTiBreakendCount > 2)
        {
            // cannot just be a DSB with TIs, so early exit for further analysis
            mType = ARM_CL_COMPLEX_OTHER;
            return;
        }
        */

        List<SvBreakend> unlinkedBreakends = mBreakends.stream().filter(x -> !tiBreakends.contains(x)).collect(Collectors.toList());

        // check for DSBs which were masked by TIs in between
        for (int i = 0; i < unlinkedBreakends.size() - 1; ++i)
        {
            final SvBreakend be1 = unlinkedBreakends.get(i);
            final SvBreakend be2 = unlinkedBreakends.get(i+1);

            if(be1.orientation() == 1 && be2.orientation() == -1)
            {
                ++dsbCount;
            }
            else if(be1.orientation() == -1 && be2.orientation() == 1
            && be2.position() - be1.position() < getMinTemplatedInsertionLength(be1, be2))
            {
                ++dsbCount;
            }
        }

        if(nonTiBreakendCount == 1)
        {
            // if after removing all TIs, all that is left is a single breakend then classify it as such
            mType = ARM_CL_ISOLATED_BE;
        }
        else if(foldbackCount == 1 && dsbCount == 0)
        {
            mType = ARM_CL_FOLDBACK;
        }
        else if(nonTiBreakendCount == 3 && foldbackCount == 1 && dsbCount == 1)
        {
            mType = ARM_CL_FOLDBACK_DSB;
        }
        else if(foldbackCount == 0 && dsbCount == 1 && nonTiBreakendCount == 2)
        {
            mType = ARM_CL_DSB;
        }
        else if(unlinkedBreakends.size() == 2 && unlinkedBreakends.get(0).orientation() == unlinkedBreakends.get(1).orientation())
        {
            mType = ARM_CL_SAME_ORIENT;
        }
        else
        {
            mType = ARM_CL_COMPLEX_OTHER;
        }
    }

    public static int[] getArmClusterData(final SvCluster cluster)
    {
        int[] results = new int[ARM_CL_MAX+1];

        for(final SvArmCluster armCluster : cluster.getArmClusters())
        {
            ++results[armCluster.getType()];
        }

        return results;
    }

    public static void logArmClusterData(final SvCluster cluster)
    {
        for(final SvArmCluster armCluster : cluster.getArmClusters())
        {
            LOGGER.debug("cluster({}) armCluster({}) breakends({}) type({})",
                    cluster.id(), armCluster.toString(), armCluster.getBreakends().size(), typeToString(armCluster.getType()));
        }
    }
}
