package com.hartwig.hmftools.isofox.novel.cohort;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_PAIR;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionFile;

public class AltSjCohortData
{
    public final AltSpliceJunctionFile AltSJ;

    // cohort data
    private final Set<String> mSampleIds;

    private Map<String,List<String>> mCancerSampleIds;

    private int mTotalFragmentCount;

    private int mMaxFragmentCount;
    private final int[] mPositionCounts; // counts at the start and end

    public AltSjCohortData(final AltSpliceJunctionFile altSJ)
    {
        AltSJ = altSJ;
        mSampleIds = Sets.newHashSet();
        mCancerSampleIds = null;

        mTotalFragmentCount = 0;

        mMaxFragmentCount = 0;
        mPositionCounts = new int[SE_PAIR];
    }

    public void addSampleCount(final String sampleId, int fragCount, final String cancerType)
    {
        if(mSampleIds.contains(sampleId))
            return;

        mSampleIds.add(sampleId);
        mTotalFragmentCount += fragCount;
        mMaxFragmentCount = max(mMaxFragmentCount, fragCount);

        if(cancerType != null)
        {
            if(mCancerSampleIds == null)
                mCancerSampleIds = Maps.newHashMap();

            List<String> samples = mCancerSampleIds.get(cancerType);

            if(samples == null)
                mCancerSampleIds.put(cancerType, Lists.newArrayList(sampleId));
            else
                samples.add(sampleId);
        }
    }

    public final Set<String> getSampleIds() { return mSampleIds; }
    public final Map<String,List<String>> cancerSampleIds() { return mCancerSampleIds; }

    public int getMaxFragmentCount() { return mMaxFragmentCount; }

    public double getAvgFragmentCount() { return mSampleIds.size() > 0 ? mTotalFragmentCount / (double)mSampleIds.size() : 0; }

    public int getPositionCount(int seIndex) { return mPositionCounts[seIndex]; }
    public void addPositionCount(int seIndex, int count) { mPositionCounts[seIndex] += count; }

    public String sampleIdsStr()
    {
        StringJoiner sj = new StringJoiner(ITEM_DELIM);
        mSampleIds.forEach(x -> sj.add(x));
        return sj.toString();
    }


}
