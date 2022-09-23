package com.hartwig.hmftools.bammetrics;

import static java.lang.Math.min;
import static java.lang.String.format;

import java.util.StringJoiner;

public class Metrics
{
    public final long[] FilterTypeCounts;
    public final int[] CoverageFrequency;

    private long mTotalFiltered;
    private Statistics mStatistics;

    public Metrics(int maxCoverage)
    {
        FilterTypeCounts = new long[FilterType.values().length];
        CoverageFrequency = new int[maxCoverage + 1];
        mTotalFiltered = -1;
        mStatistics = null;
    }

    public long totalFiltered() { return mTotalFiltered; }
    public Statistics statistics() { return mStatistics;}

    public double calcFilteredPercentage(final FilterType type)
    {
        calcTotalFiltered();

        if(mTotalFiltered == 0)
            return 0;

        return FilterTypeCounts[type.ordinal()] / (double)mTotalFiltered;
    }

    public double calcCoverageFrequency(int coverageLevel)
    {
        calcTotalFiltered();

        if(mTotalFiltered == 0)
            return 0;

        int frequencyTotal = 0;

        for(int i = 0; i < min(CoverageFrequency.length, coverageLevel); ++i)
        {
            frequencyTotal += CoverageFrequency[i];
        }

        return frequencyTotal / (double)mTotalFiltered;
    }

    private void calcStatistics()
    {
        if(mStatistics != null)
            return;



        mStatistics = new Statistics(0, 0, 0);
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
