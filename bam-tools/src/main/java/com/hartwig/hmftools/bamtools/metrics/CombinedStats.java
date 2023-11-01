package com.hartwig.hmftools.bamtools.metrics;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;

public class CombinedStats
{
    private final CoverageMetrics mCoverageMetrics;
    private final FragmentLengths mFragmentLengths;
    private final ReadCounts mReadCounts;
    private final FlagStats mFlagStats;

    private final List<TargetRegionStats> mTargetRegions;

    private final PerformanceCounter mPerfCounter;

    public CombinedStats(int maxCoverage)
    {
        mCoverageMetrics = new CoverageMetrics(maxCoverage);
        mFragmentLengths = new FragmentLengths();
        mPerfCounter = new PerformanceCounter("Coverage");
        mReadCounts = new ReadCounts();
        mFlagStats = new FlagStats();
        mTargetRegions = Lists.newArrayList();
    }

    public CoverageMetrics coverageMetrics() { return mCoverageMetrics; }
    public FragmentLengths fragmentLengths() { return mFragmentLengths; }

    public PerformanceCounter perfCounter() { return mPerfCounter; }
    public ReadCounts readCounts() { return mReadCounts; }

    public FlagStats flagStats()
    {
        return mFlagStats;
    }
    public List<TargetRegionStats> targetRegions() { return mTargetRegions; }

    public synchronized void addStats(
            final CoverageMetrics metrics, final FragmentLengths fragmentLengths, final ReadCounts readCounts,
            final FlagStats flagStats, final List<TargetRegionStats> targetRegions, final PerformanceCounter perfCounter)
    {
        mCoverageMetrics.merge(metrics);
        mFragmentLengths.merge(fragmentLengths);
        mReadCounts.merge(readCounts);
        mPerfCounter.merge(perfCounter);
        mFlagStats.merge(flagStats);

        for(TargetRegionStats targetRegion : targetRegions)
        {
            TargetRegionStats matchedRegion = mTargetRegions.stream().filter(x -> x.Region.matches(targetRegion.Region)).findFirst().orElse(null);

            if(matchedRegion != null)
            {
                matchedRegion.Counts.merge(targetRegion.Counts);
            }
            else
            {
                mTargetRegions.add(targetRegion);
            }
        }
    }
}
