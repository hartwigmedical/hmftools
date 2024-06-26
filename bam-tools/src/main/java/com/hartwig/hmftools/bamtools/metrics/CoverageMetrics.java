package com.hartwig.hmftools.bamtools.metrics;

import static java.lang.Math.abs;
import static java.lang.Math.pow;
import static java.lang.Math.round;
import static java.lang.Math.sqrt;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;

import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

public class CoverageMetrics
{
    public final long[] FilterTypeCounts;
    public final long[] CoverageFrequency;

    private long mCoverageBases; // bases with any level of coverage

    private long mTotalBases; // filtered and unfiltered
    private Statistics mStatistics;

    public CoverageMetrics(int maxCoverage)
    {
        FilterTypeCounts = new long[FilterType.values().length];
        CoverageFrequency = new long[maxCoverage + 1];
        mCoverageBases = 0;
        mTotalBases = -1;
        mStatistics = null;
    }

    public void addCoverageBases(long bases) { mCoverageBases += bases; }
    public long coverageBases() { return mCoverageBases; }
    public long zeroCoverageBases() { return CoverageFrequency[0]; }

    // note this does not include unmappable regions (eg 'N' in centromeres and telomeres)
    public long genomeTerritory() { return zeroCoverageBases() + coverageBases(); }

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

    public double calcCoverageFrequency(int coverageLevel, double genomeTerritory)
    {
        if(genomeTerritory <= 0)
            return 0;

        // calculates the percentage of bases with this coverage level or higher, exclude the bases with zero coverage

        long frequencyTotal = 0;

        for(int i = 0; i < CoverageFrequency.length; ++i)
        {
            if(i < coverageLevel)
                continue;

            frequencyTotal += CoverageFrequency[i];
        }

        return frequencyTotal / genomeTerritory;
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

            varianceTotal += pow(i - mean, 2) * CoverageFrequency[i];
        }

        double stdDeviation = sqrt(varianceTotal / total);

        // calculate median difference from median
        List<DiffFrequency> diffFrequencies = Lists.newArrayList();
        for(int i = 0; i < CoverageFrequency.length; ++i)
        {
            int coverageDiff = (int)round(abs(i - median));
            addMedianDiffFrequency(diffFrequencies, coverageDiff, CoverageFrequency[i]);
        }

        cumulativeTotal = 0;
        int mad = 0;
        for(DiffFrequency diffFrequency : diffFrequencies)
        {
            if(cumulativeTotal + diffFrequency.Frequency >= medianTotal)
            {
                mad = diffFrequency.Diff;
                break;
            }
            else
            {
                cumulativeTotal += diffFrequency.Frequency;
            }
        }

        BT_LOGGER.debug(format("mean(%.3f) median(%.1f) mad(%d) totalFrequency=%d total=%d) stdDeviation(%.3f)",
                mean, median, mad, totalFrequency, total, stdDeviation));

        mStatistics = new Statistics(mean, median, mad, stdDeviation);
    }

    private class DiffFrequency
    {
        public final int Diff;
        public long Frequency;

        public DiffFrequency(int diff, long frequency)
        {
            Diff = diff;
            Frequency = frequency;
        }
    }

    private void addMedianDiffFrequency(final List<DiffFrequency> diffFrequencies, int diff, long frequency)
    {
        int index = 0;
        while(index < diffFrequencies.size())
        {
            if(diffFrequencies.get(index).Diff == diff)
            {
                diffFrequencies.get(index).Frequency += frequency;
                return;
            }
            else if(diffFrequencies.get(index).Diff > diff)
            {
                break;
            }

            ++index;
        }

        diffFrequencies.add(index, new DiffFrequency(diff, frequency));
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

    public synchronized void merge(final CoverageMetrics other)
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

    public static int getCoverageBucket(int coverage)
    {
        // round to nearest unit up to 1000, then 10s up to 3000 then 100s
        if(coverage <= 100)
            return coverage;

        if(coverage <= 1000)
            return 10 * (int)round(coverage/10.0);

        if(coverage <= 10000)
            return 100 * (int)round(coverage/100.0);

        return 1000;
    }
}
