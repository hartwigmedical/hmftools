package com.hartwig.hmftools.svanalysis.types;

import static java.lang.Math.max;

import java.util.List;

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
    private final String mArm;
    private long mStartPos;
    private long mEndPos;
    private double mMinCopyNumber;
    private double mMaxCopyNumber;

    private int mType;
    private int mTICount;

    private static final Logger LOGGER = LogManager.getLogger(SvArmCluster.class);

    public SvArmCluster(int id, final SvCluster cluster, final String chr, final String arm)
    {
        mId = id;
        mCluster = cluster;
        mArm = arm;
        mChromosome = chr;
        mStartPos = 0;
        mEndPos = 0;
        mBreakends = Lists.newArrayList();
        mMaxCopyNumber = 0;
        mMinCopyNumber = 0;
        mType = -1;
        mTICount = 0;
    }

    public int id() { return mId; }
    public final String toString()
    {
        return String.format("%d: %s %s_%s %d:%d",
                mId, typeToString(mType), mChromosome, mArm, mStartPos, mEndPos);
    }

    public final String chromosome() { return mChromosome; }
    public final String arm() { return mArm; }
    public long posStart() { return mStartPos; }
    public long posEnd() { return mEndPos; }

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

        double lowCopyNumber = breakend.getCopyNumber(true);

        if(mMinCopyNumber == 0 || lowCopyNumber < mMinCopyNumber)
            mMinCopyNumber = lowCopyNumber;

        mMaxCopyNumber = max(mMaxCopyNumber, breakend.getCopyNumber(false));
    }

    public double getMinCopyNumber() { return mMinCopyNumber; }
    public double getMaxCopyNumber() { return mMaxCopyNumber; }

    public static final int ARM_CL_ISOLATED_BE = 0;
    public static final int ARM_CL_TI_ONLY = 1;
    public static final int ARM_CL_DSB = 2;
    public static final int ARM_CL_FOLDBACK = 3;
    public static final int ARM_CL_FOLDBACK_DSB = 4;
    public static final int ARM_CL_COMPLEX_FOLDBACK = 5;
    public static final int ARM_CL_COMPLEX_LINE = 6;
    public static final int ARM_CL_COMPLEX_OTHER = 7;

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
        }

        return "Unknown";
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

            if(var1.getFoldbackLink(be1.usesStart()).equals(var2.id()))
            {
                mType = ARM_CL_FOLDBACK;
                return;
            }

            if(var1 == var2)
            {
                mType = ARM_CL_ISOLATED_BE;
                return;
            }

            final SvLinkedPair tiPair = var1.getLinkedPair(be1.usesStart());

            if(tiPair != null && tiPair.hasBreakend(be2, true))
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

        // otherwise count the number of foldbacks, TIs, DSBs and consecutive BEs to determine the type
        int dsbCount = 0;
        int foldbackCount = 0;
        int consecCount = 0;
        int suspectLine = 0;
        int tiCount = 0;

        for (int i = 0; i < mBreakends.size() - 1; ++i)
        {
            final SvBreakend be1 = mBreakends.get(i);
            final SvBreakend be2 = mBreakends.get(i+1);

            if(be1.getSV().isLineElement(be1.usesStart()))
            {
                ++suspectLine;
            }

            final SvLinkedPair dbPair = be1.getSV().getDBLink(be1.usesStart());

            if(dbPair != null && dbPair == be2.getSV().getDBLink(be2.usesStart()))
            {
                ++dsbCount;
            }

            final SvLinkedPair tiPair = be1.getSV().getLinkedPair(be1.usesStart());

            if(tiPair != null && tiPair.hasBreakend(be2, true))
            {
                ++tiCount;
            }

            if(be1.getSV().getFoldbackLink(be1.usesStart()).equals(be2.getSV().id()))
            {
                ++foldbackCount;
            }
            else if(be1.getSV().getConsecBEStart(be1.usesStart()).equals(be2.getSV().id()))
            {
                ++consecCount;
            }
        }

        mTICount = tiCount;

        if(mBreakends.size() == tiCount * 2 + 1)
        {
            // if after removing all TIs, all that is left is a single breakend then classify it as such
            mType = ARM_CL_ISOLATED_BE;
        }

        if(suspectLine > 0)
        {
            mType = ARM_CL_COMPLEX_LINE;
        }
        else if(mBreakends.size() == 3 && foldbackCount == 1 && dsbCount == 1)
        {
            mType = ARM_CL_FOLDBACK_DSB;
        }
        else if(foldbackCount == 0 && consecCount == 0 && dsbCount >= mBreakends.size() / 2)
        {
            mType = ARM_CL_DSB;
        }
        else if(foldbackCount > 0)
        {
            mType = ARM_CL_COMPLEX_FOLDBACK;
        }
        else
        {
            mType = ARM_CL_COMPLEX_OTHER;
        }
    }

    public static int[] getArmClusterData(final SvCluster cluster)
    {
        int[] results = new int[ARM_CL_COMPLEX_OTHER+1];

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
