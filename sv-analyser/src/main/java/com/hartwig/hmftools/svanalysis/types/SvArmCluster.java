package com.hartwig.hmftools.svanalysis.types;

import static java.lang.Math.max;

import static com.hartwig.hmftools.svanalysis.types.SvCluster.isSpecificCluster;

import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.LoggerConfig;

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
    }

    public int id() { return mId; }
    public final String toString()
    {
        return String.format("%d: %s_%s %d:%d", mId, mChromosome, mArm, mStartPos, mEndPos);
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

    public static final int ARM_CL_SINGLE = 0;
    public static final int ARM_CL_REMOTE_TI = 1;
    public static final int ARM_CL_DSB = 2;
    public static final int ARM_CL_MULTIPLE_DSBS = 3;
    public static final int ARM_CL_FOLDBACK = 4;
    public static final int ARM_CL_FOLDBACK_DSB = 5;
    public static final int ARM_CL_FOLDBACK_TI = 6;
    public static final int ARM_CL_FOLDBACK_PAIR_FACING = 7;
    public static final int ARM_CL_FOLDBACK_PAIR_OPPOSING = 8;
    public static final int ARM_CL_FOLDBACK_PAIR_SAME = 9;
    public static final int ARM_CL_COMPLEX_FOLDBACK = 10;
    public static final int ARM_CL_COMPLEX_LINE = 11;
    public static final int ARM_CL_COMPLEX_OTHER = 12;

    public static final String typeToString(int type)
    {
        switch(type)
        {
            case ARM_CL_SINGLE : return "SINGLE";
            case ARM_CL_REMOTE_TI : return "REMOTE_TI";
            case ARM_CL_DSB : return "DSB";
            case ARM_CL_MULTIPLE_DSBS : return "MULTIPLE_DSB";
            case ARM_CL_FOLDBACK : return "FOLDBACK";
            case ARM_CL_FOLDBACK_DSB : return "FOLDBACK_DSB";
            case ARM_CL_FOLDBACK_TI : return "FOLDBACK_TI";
            case ARM_CL_FOLDBACK_PAIR_SAME : return "FOLDBACK_PAIR_SAME";
            case ARM_CL_FOLDBACK_PAIR_FACING : return "FOLDBACK_PAIR_OPP";
            case ARM_CL_FOLDBACK_PAIR_OPPOSING : return "FOLDBACK_PAIR_FACE";
            case ARM_CL_COMPLEX_FOLDBACK : return "COMPLEX_FOLDBACK";
            case ARM_CL_COMPLEX_LINE: return "COMPLEX_LINE";
            case ARM_CL_COMPLEX_OTHER : return "COMPLEX_OTHER";
        }

        return "Unknown";
    }

    public int getType()
    {
        if(mBreakends.size() == 1)
            return ARM_CL_SINGLE;

        if(mBreakends.size() == 2)
        {
            final SvBreakend be1 = mBreakends.get(0);
            final SvVarData var1 = be1.getSV();
            final SvBreakend be2 = mBreakends.get(1);
            final SvVarData var2 = be2.getSV();

            if(var1.getFoldbackLink(be1.usesStart()).equals(var2.id()))
                return ARM_CL_FOLDBACK;

            if(var1 == var2)
                return ARM_CL_SINGLE;

            final SvLinkedPair tiPair = var1.getLinkedPair(be1.usesStart());

            if(tiPair != null && tiPair == var2.getLinkedPair(be2.usesStart()))
            {
                return ARM_CL_REMOTE_TI;
            }
        }

        // otherwise count the number of foldbacks, DSBs and consecutive BEs to determine the type
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

            if(tiPair != null && tiPair == be2.getSV().getLinkedPair(be2.usesStart()))
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

        if(suspectLine > 0)
            return ARM_CL_COMPLEX_LINE;

        // some special sub-groups
        if(mBreakends.size() == 3 && foldbackCount == 1)
        {
            if(dsbCount == 1)
                return ARM_CL_FOLDBACK_DSB;

            if(tiCount == 1)
                return ARM_CL_FOLDBACK_TI;
        }
        else if(mBreakends.size() == 4 && foldbackCount == 2)
        {
            final SvBreakend be1 = mBreakends.get(0);
            final SvBreakend be4 = mBreakends.get(3);

            if(be1.orientation() == be4.orientation())
            {
                return ARM_CL_FOLDBACK_PAIR_SAME;
            }
            else if(be1.orientation() == -1 && be4.orientation() == 1)
            {
                return ARM_CL_FOLDBACK_PAIR_FACING;
            }
            else
            {
                return ARM_CL_FOLDBACK_PAIR_OPPOSING;
            }
        }

        if(foldbackCount == 0 && consecCount == 0 && dsbCount >= mBreakends.size() / 2)
            return mBreakends.size() == 2 ? ARM_CL_DSB : ARM_CL_MULTIPLE_DSBS;

        if(foldbackCount > 0)
            return ARM_CL_COMPLEX_FOLDBACK;

        return ARM_CL_COMPLEX_OTHER;
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

    public static void mergeArmClusters(List<SvArmCluster> armClusters)
    {
        // merge if any 2 are linked by the same foldback
        for(int i = 0; i < armClusters.size(); ++i)
        {
            SvArmCluster ac1 = armClusters.get(i);
            List<SvBreakend> bl1 = ac1.getBreakends();

            int j = i+1;
            while(j < armClusters.size())
            {
                SvArmCluster ac2 = armClusters.get(j);

                if(!ac1.chromosome().equals(ac2.chromosome()) || ac1.arm() != ac2.arm())
                {
                    ++j;
                    continue;
                }

                // check for a matching foldback
                List<SvBreakend> bl2 = ac2.getBreakends();

                boolean foldbackFound = false;

                for(final SvBreakend be1 : bl1)
                {
                    for(final SvBreakend be2 : bl2)
                    {
                        if(be1.getSV().getFoldbackLink(be1.usesStart()).equals(be2.getSV().id()))
                        {
                            foldbackFound = true;
                            break;
                        }
                    }

                    if(foldbackFound)
                        break;
                }

                if(foldbackFound)
                {
                    for(final SvBreakend be2 : bl2)
                    {
                        ac1.addBreakend(be2);
                    }

                    armClusters.remove(j);
                }
                else
                {
                    ++j;
                }
            }
        }
    }
}
