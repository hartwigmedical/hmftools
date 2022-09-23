package com.hartwig.hmftools.bammetrics;

import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import static java.lang.String.format;

import java.util.StringJoiner;

public class Metrics
{
    public final long[] FilterTypeCounts;
    public final int[] CoverageFrequency;

    private int mCoverageBases; // bases with any level of coverage

    private long mTotalFiltered;
    private Statistics mStatistics;

    public Metrics(int maxCoverage)
    {
        FilterTypeCounts = new long[FilterType.values().length];
        CoverageFrequency = new int[maxCoverage + 1];
        mCoverageBases = 0;
        mTotalFiltered = -1;
        mStatistics = null;
    }

    public void addCoverageBases(int bases) { mCoverageBases += bases; }
    public int coverageBases() { return mCoverageBases; }
    public int zeroCoverageBases() { return CoverageFrequency[0]; }

    public Statistics statistics() { return mStatistics; }

    public void finalise(boolean excludeZeroCoverge)
    {
        calcTotalFiltered();

        calcStatistics(excludeZeroCoverge);
    }

    public double calcFilteredPercentage(final FilterType type)
    {
        calcTotalFiltered();

        if(mTotalFiltered == 0)
            return 0;

        return FilterTypeCounts[type.ordinal()] / (double)mTotalFiltered;
    }

    public double calcCoverageFrequency(int coverageLevel)
    {
        long totalCoverage = FilterTypeCounts[FilterType.UNFILTERED.ordinal()];

        long frequencyTotal = 0;

        for(int i = 0; i < CoverageFrequency.length; ++i)
        {
            if(i < coverageLevel)
                continue;

            frequencyTotal += CoverageFrequency[i] * i;
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
            totalFrequency += CoverageFrequency[i] * i;
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

            varianceTotal += CoverageFrequency[i] * pow(i - mean, 2);
        }

        double stdDeviation = sqrt(varianceTotal / total);
        mStatistics = new Statistics(mean, median, stdDeviation, 0);
    }

    private void calcTotalFiltered()
    {
        if(mTotalFiltered >= 0)
            return;

        mTotalFiltered = 0;

        for(FilterType type : FilterType.values())
        {
            if(type == FilterType.UNFILTERED)
                continue;

            mTotalFiltered += FilterTypeCounts[type.ordinal()];
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
