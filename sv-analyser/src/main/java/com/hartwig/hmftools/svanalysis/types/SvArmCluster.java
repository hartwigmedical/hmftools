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
    public static int ARM_CL_SIMPLE_FOLDBACK = 4;
    public static int ARM_CL_COMPLEX_FOLDBACK = 5;
    public static int ARM_CL_COMPLEX_OTHER = 6;

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
                return ARM_CL_SIMPLE_FOLDBACK;

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

        // otherwise count the number of foldbacks, DSBs and consecutive BEs to determine the type
        int dsbCount = 0;
        int foldbackCount = 0;
        int consecCount = 0;

        for (int i = 0; i < mBreakends.size() - 1; ++i)
        {
            final SvBreakend be1 = mBreakends.get(i);
            final SvBreakend be2 = mBreakends.get(i+1);

            final SvLinkedPair dbPair = be1.getSV().getDBLink(be1.usesStart());

            if(dbPair != null && dbPair == be2.getSV().getDBLink(be2.usesStart()))
            {
                ++dsbCount;
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

        if(foldbackCount == 0 && consecCount == 0 && dsbCount >= mBreakends.size() / 2)
            return mBreakends.size() == 2 ? ARM_CL_DSB : ARM_CL_MULTIPLE_DSBS;

        if(foldbackCount > 0)
            return ARM_CL_COMPLEX_FOLDBACK;

        return ARM_CL_COMPLEX_OTHER;
    }

    public static int[] getArmClusterData(final SvCluster cluster)
    {
        // isSpecificCluster(cluster);

        int[] results = new int[ARM_CL_COMPLEX_OTHER+1];

        for(final SvArmCluster armCluster : cluster.getArmClusters())
        {
            ++results[armCluster.getType()];
        }

        return results;
    }

    public static void mergeArmClusters(List<SvArmCluster> armClusters)
    {
        // merge if any 2 have SV in a foldback
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
