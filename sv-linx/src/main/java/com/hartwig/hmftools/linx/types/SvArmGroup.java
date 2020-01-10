package com.hartwig.hmftools.linx.types;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.calcConsistency;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.makeChrArmStr;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;

import com.google.common.collect.Lists;

import java.util.List;

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

    public void populateDbData(final ClusterMetrics data)
    {
        final List<SvLinkedPair> processedDBs = Lists.newArrayList();

        for(int i = 0; i < mBreakends.size(); ++i)
        {
            final SvBreakend breakend = mBreakends.get(i);

            if(breakend.getSV().type() == DEL && breakend.orientation() == 1)
            {
                // check for stand-alone DELs
                if(i < mBreakends.size() - 1 && mBreakends.get(i + 1).getSV() == breakend.getSV())
                {
                    ++data.DBCount;
                    ++data.ClusterDBCount;

                    data.TotalDBLength += breakend.getSV().length();

                    if(breakend.getSV().length() < 100)
                        ++data.ShortDBCount;

                    ++i;
                    continue;
                }
            }

            final SvLinkedPair dbLink = breakend.getDBLink();

            if(dbLink == null || processedDBs.contains(dbLink))
                continue;

            processedDBs.add(dbLink);
            ++data.DBCount;

            if(dbLink.length() < 100)
                ++data.ShortDBCount;

            if(dbLink.getOtherBreakend(breakend).getCluster() == mCluster)
            {
                ++data.ClusterDBCount;
                data.TotalDBLength += max(dbLink.length(), 0);
            }
        }
    }

}
