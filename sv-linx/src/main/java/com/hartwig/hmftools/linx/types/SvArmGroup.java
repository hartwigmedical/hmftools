package com.hartwig.hmftools.linx.types;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.linx.analysis.SvUtilities.calcConsistency;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.makeChrArmStr;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.isStart;

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

    // width of SVs on the arm taking into account any excluded SVs
    private long mStartPos;
    private long mEndPos;

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

    public boolean matches(final SvArmGroup other)
    {
        return mChromosome.equals(other.chromosome()) && mArm.equals(other.arm());
    }

    public int consistency()
    {
        return mBreakends.stream().mapToInt(x -> calcConsistency(x)).sum();
    }

    public boolean isConsistent() { return consistency() == 0; }

    public static final int DB_DATA_COUNT = 0;
    public static final int DB_DATA_CLUSTER_COUNT = 1;
    public static final int DB_DATA_SHORT_COUNT = 2;
    public static final int DB_DATA_TOTAL_LENGTH = 3;
    public static final int DB_DATA_BOUNDARY_LENGTH = 4;

    public void populateDbData(final int[] data)
    {
        List<SvLinkedPair> processedDBs = Lists.newArrayList();

        for(final SvBreakend breakend : mBreakends)
        {
            final SvLinkedPair dbLink = breakend.getDBLink();

            if(dbLink == null || processedDBs.contains(dbLink))
                continue;

            processedDBs.add(dbLink);
            ++data[DB_DATA_COUNT];

            if(dbLink.length() < 100)
                ++data[DB_DATA_SHORT_COUNT];

            if(dbLink.getOtherBreakend(breakend).getCluster() == mCluster)
            {
                ++data[DB_DATA_CLUSTER_COUNT];
                data[DB_DATA_TOTAL_LENGTH] += max(dbLink.length(), 0);
            }
        }

        // if this cluster has DBs between variants on the arm, then report the total width covered by the breakends
        if(data[DB_DATA_CLUSTER_COUNT] > 0)
        {
            data[DB_DATA_BOUNDARY_LENGTH] += (int)(mEndPos - mStartPos);
        }
    }

}
