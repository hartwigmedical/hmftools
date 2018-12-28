package com.hartwig.hmftools.svanalysis.types;

import static com.hartwig.hmftools.svanalysis.types.SvCluster.isSpecificCluster;

import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SvArmCluster
{
    private int mId;
    private List<SvBreakend> mBreakends;
    final SvCluster mCluster;

    private final String mChromosome;
    private final String mArm;
    private long mStartPos;
    private long mEndPos;

    private static final Logger LOGGER = LogManager.getLogger(SvArmGroup.class);

    public SvArmCluster(int id, final SvCluster cluster, final String chr, final String arm)
    {
        mCluster = cluster;
        mArm = arm;
        mChromosome = chr;
        mStartPos = 0;
        mEndPos = 0;
        mBreakends = Lists.newArrayList();
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
    }

    public static int ARM_CL_SINGLE = 0;
    public static int ARM_CL_REMOTE_TI = 1;
    public static int ARM_CL_DSB = 2;
    public static int ARM_CL_MULTIPLE_DSBS = 3;
    public static int ARM_CL_FOLDBACK = 4;
    public static int ARM_CL_COMPLEX = 5;

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
            {
                return ARM_CL_SINGLE;
            }
            else
            {
                final SvLinkedPair tiPair = var1.getLinkedPair(be1.usesStart());

                if(tiPair != null && tiPair == var2.getLinkedPair(be2.usesStart()))
                {
                    return ARM_CL_REMOTE_TI;
                }
            }
        }

        // look for all breakends forming DSBs with each other
        if((mBreakends.size() % 2) == 0)
        {
            boolean allDsbs = true;

            for (int i = 0; i < mBreakends.size() - 1; i = i+2)
            {
                final SvBreakend be1 = mBreakends.get(i);
                final SvBreakend be2 = mBreakends.get(i+1);

                final SvLinkedPair dbPair = be1.getSV().getDBLink(be1.usesStart());

                if(dbPair == null || dbPair != be2.getSV().getDBLink(be2.usesStart()))
                {
                    allDsbs = false;
                }
            }

            if(allDsbs)
                return mBreakends.size() == 2 ? ARM_CL_DSB : ARM_CL_MULTIPLE_DSBS;
        }

        return ARM_CL_COMPLEX;
    }

    public static int[] getArmClusterData(final SvCluster cluster)
    {
        // isSpecificCluster(cluster);

        int[] results = new int[ARM_CL_COMPLEX+1];

        for(final SvArmCluster armCluster : cluster.getArmClusters())
        {
            ++results[armCluster.getType()];
        }

        return results;
    }
}
