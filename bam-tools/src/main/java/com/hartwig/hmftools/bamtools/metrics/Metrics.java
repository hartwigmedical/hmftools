package com.hartwig.hmftools.bamtools.metrics;

import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;

import java.util.StringJoiner;

public class Metrics
{
    public final long[] FilterTypeCounts;
    public final int[] CoverageFrequency;

    private long mCoverageBases; // bases with any level of coverage

    private long mTotalBases; // filtered and unfiltered
    private Statistics mStatistics;

    public Metrics(int maxCoverage)
    {
        FilterTypeCounts = new long[FilterType.values().length];
        CoverageFrequency = new int[maxCoverage + 1];
        mCoverageBases = 0;
        mTotalBases = -1;
        mStatistics = null;
    }

    public void addCoverageBases(long bases) { mCoverageBases += bases; }
    public long coverageBases() { return mCoverageBases; }
    public long zeroCoverageBases() { return CoverageFrequency[0]; }

    public Statistics statistics() { return mStatistics; }

    public void finalise(boolean excludeZeroCoverge)
    {
        calcTotalFiltered();

        calcStatistics(excludeZeroCoverge);
    }

    public double calcFilteredPercentage(final FilterType type)
    {
        // calculate percentage of a filtered type vs all bases including unfiltered
        calcTotalFiltered();

        if(mTotalBases == 0)
            return 0;

        return FilterTypeCounts[type.ordinal()] / (double) mTotalBases;
    }

    public double calcCoverageFrequency(int coverageLevel)
    {
        // calculates the percentage of bases with this coverage level or higher, exclude the bases with zero coverage
        long totalCoverage = FilterTypeCounts[FilterType.UNFILTERED.ordinal()];

        long frequencyTotal = 0;

        for(int i = 0; i < CoverageFrequency.length; ++i)
        {
            if(i < coverageLevel)
                continue;

            frequencyTotal += (long)CoverageFrequency[i] * i;
        }

        return frequencyTotal / (double)totalCoverage;
    }

    private void calcStatistics(boolean excludeZeroCoverge)
    {
        if(mStatistics != null)
            return;

        long total = 0;
        long totalFrequency = 0;
        for(int i = 0; i < CoverageFrequency.length; ++i)
        {
            if(excludeZeroCoverge && i == 0)
                continue;

            total += CoverageFrequency[i];
            totalFrequency += (long)CoverageFrequency[i] * i;
        }

        if(total == 0)
        {
            mStatistics = new Statistics(0, 0, 0, 0);
            return;
        }

        double mean = totalFrequency / (double)total;

        double median = -1;
        double varianceTotal = 0;
        long cumulativeTotal = 0;
        long medianTotal = total / 2;
        for(int i = 0; i < CoverageFrequency.length; ++i)
        {
            if(excludeZeroCoverge && i == 0)
                continue;

            if(median < 0)
            {
                if(cumulativeTotal + CoverageFrequency[i] >= medianTotal)
                    median = i;
                else
                    cumulativeTotal += CoverageFrequency[i];
            }

            varianceTotal += pow(i - mean, 2) * CoverageFrequency[i];
        }

        double stdDeviation = sqrt(varianceTotal / total);

        BM_LOGGER.debug(format("mean(%.3f totalFrequency=%d total=%d) stdDeviation(%.3f variantTotal=%.3f)",
                mean, totalFrequency, total, stdDeviation, varianceTotal));

        mStatistics = new Statistics(mean, median, stdDeviation, 0);
    }

    private void calcTotalFiltered()
    {
        if(mTotalBases >= 0)
            return;

        mTotalBases = 0;

        for(FilterType type : FilterType.values())
        {
            mTotalBases += FilterTypeCounts[type.ordinal()];
        }
    }

    public synchronized void merge(final Metrics other)
    {
        for(FilterType type : FilterType.values())
        {
            FilterTypeCounts[type.ordinal()] += other.FilterTypeCounts[type.ordinal()];
        }

        for(int i = 0; i < CoverageFrequency.length; ++i)
        {
            CoverageFrequency[i] += other.CoverageFrequency[i];
        }

        addCoverageBases(other.coverageBases());
    }

    public String toString()
    {
        StringJoiner sj = new StringJoiner(" ");
        for(FilterType type : FilterType.values())
        {
            sj.add(format("%s=%d", type, FilterTypeCounts[type.ordinal()]));
        }

        return format("counts(%s)", sj);
    }
}
