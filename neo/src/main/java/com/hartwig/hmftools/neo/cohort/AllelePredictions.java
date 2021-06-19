package com.hartwig.hmftools.neo.cohort;

import static java.lang.Math.round;

import static com.hartwig.hmftools.common.stats.Percentiles.PERCENTILE_COUNT;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;

public class AllelePredictions
{
    public final String Allele;

    private final Set<String> mPeptides; // to avoid double-counting

    private final int[] mFrequencies;
    private int mTotal;

    private static final int AFFINITY_LVL_1 = 1000; //
    private static final int AFFINITY_MAX = 50000;
    private static final int AFFINITY_STEP = 10;

    public AllelePredictions(final String allele)
    {
        Allele = allele;
        mPeptides = Sets.newHashSet();
        mFrequencies = new int[calcMaxBucket()];
        mTotal = 0;
    }

    public int[] getFrequencies() { return mFrequencies; }
    public int getTotal() { return mTotal; }

    private int calcMaxBucket()
    {
        return AFFINITY_LVL_1 + (AFFINITY_MAX - AFFINITY_LVL_1) / AFFINITY_STEP;
    }

    private int bucket(double affinity)
    {
        if(affinity <= AFFINITY_LVL_1)
            return (int)round(affinity);

        if(affinity >= AFFINITY_MAX)
            return mFrequencies.length - 1;

        return AFFINITY_LVL_1 + (int)round((affinity - AFFINITY_LVL_1) / 10);
    }

    private int affinity(int bucket)
    {
        if(bucket <= AFFINITY_LVL_1)
            return bucket;

        if(bucket >= mFrequencies.length)
            return AFFINITY_MAX;

        return AFFINITY_LVL_1 + (bucket - AFFINITY_LVL_1) * 10;
    }

    public void addPeptide(final String peptide, double affinity)
    {
        if(mPeptides.contains(peptide))
            return;

        mPeptides.add(peptide);

        int bucket = bucket(affinity);

        if(bucket < mFrequencies.length)
            ++mFrequencies[bucket];
        else
            ++mFrequencies[mFrequencies.length - 1];

        ++mTotal;
    }

    public int[] buildPercentiles()
    {
        int[] percentileBuckets = buildPercentiles(mFrequencies, mTotal, PERCENTILE_COUNT);
        int[] percentiles = new int[PERCENTILE_COUNT];
        for(int i = 0; i < percentiles.length; ++i)
        {
            percentiles[i] = affinity(percentileBuckets[i]);
        }

        return percentiles;
    }

    public static int[] buildPercentiles(final int[] frequences, int total, int percentileCount)
    {
        int itemsPerPercentile = (int)round(total / (double)percentileCount);

        int[] percentiles = new int[percentileCount];

        int itemCount = 0;
        int percentile = 0;

        for(int i = 0; i < frequences.length; ++i)
        {
            int next = frequences[i];

            if(itemCount + next < itemsPerPercentile)
            {
                itemCount += next;
                continue;
            }

            int used = itemsPerPercentile - itemCount;

            int remainder = next - used;
            percentiles[percentile++] = i;

            while(remainder > itemsPerPercentile)
            {
                percentiles[percentile++] = i;
                remainder -= itemsPerPercentile;
            }

            itemCount = remainder;
        }

        return percentiles;
    }
}
