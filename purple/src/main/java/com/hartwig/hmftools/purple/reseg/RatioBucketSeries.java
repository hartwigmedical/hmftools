package com.hartwig.hmftools.purple.reseg;

import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_RATIO_BUCKET_MAX;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_RATIO_BUCKET_MIN;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_RATIO_BUCKET_STEP;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_ROLLING_AVG_WINDOW;

import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.Lists;

public final class RatioBucketSeries
{
    // used to cache and manipulate binned observed tumor ratios and diffs
    public record Bucket(double level, double value) {}

    private final double mRatioMin;
    private final double mRatioMax;
    private final double mRatioStep;
    private final double[] mCounts;
    private final int mBucketCount;

    public RatioBucketSeries(final double ratioMin, final double ratioMax, final double ratioStep)
    {
        mRatioMin = ratioMin;
        mRatioMax = ratioMax;
        mRatioStep = ratioStep;

        mBucketCount = (int)Math.round((ratioMax - ratioMin) / ratioStep) + 1;
        mCounts = new double[mBucketCount];
    }

    public RatioBucketSeries()
    {
        this(RESEG_RATIO_BUCKET_MIN, RESEG_RATIO_BUCKET_MAX, RESEG_RATIO_BUCKET_STEP);
    }

    public void addCount(double value, double weight)
    {
        if(value < mRatioMin || value > mRatioMax)
            return;

        int index = (int)Math.round((value - mRatioMin) / mRatioStep);
        index = Math.max(0, Math.min(mBucketCount - 1, index));
        mCounts[index] += weight;
    }

    public List<Bucket> buildSmoothedSeries() { return buildSmoothedSeries(RESEG_ROLLING_AVG_WINDOW); }

    public List<Bucket> buildSmoothedSeries(int rollingAvergeCount)
    {
        // rolling-average, edge-trimmed (buckets without a full window are dropped, then duplicates are removed
        int half = rollingAvergeCount / 2;
        List<Bucket> smoothed = Lists.newArrayList();

        for(int i = half; i < mBucketCount - half; i++)
        {
            double sum = 0;
            for(int k = i - half; k <= i + half; k++)
            {
                sum += mCounts[k];
            }

            double avg = sum / rollingAvergeCount;
            double level = Math.round((mRatioMin + i * mRatioStep) / mRatioStep) * mRatioStep;
            smoothed.add(new Bucket(level, avg));
        }

        return removeDuplicates(smoothed);
    }

    private static List<Bucket> removeDuplicates(final List<Bucket> series)
    {
        List<Bucket> collapsed = new ArrayList<>();

        for(Bucket bucket : series)
        {
            if(!collapsed.isEmpty() && collapsed.get(collapsed.size() - 1).value() == bucket.value())
                continue;

            collapsed.add(bucket);
        }

        return collapsed;
    }
}
