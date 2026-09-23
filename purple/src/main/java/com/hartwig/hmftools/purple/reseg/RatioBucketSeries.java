package com.hartwig.hmftools.purple.reseg;

import static com.hartwig.hmftools.common.utils.Doubles.round;
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

    private final double[] mCounts;
    private final int mBucketCount;

    public RatioBucketSeries()
    {
        mBucketCount = (int)Math.round((RESEG_RATIO_BUCKET_MAX - RESEG_RATIO_BUCKET_MIN) / RESEG_RATIO_BUCKET_STEP) + 1;
        mCounts = new double[mBucketCount];
    }

    public void addCount(double value, double weight)
    {
        if(value < RESEG_RATIO_BUCKET_MIN || value > RESEG_RATIO_BUCKET_MAX)
            return;

        int index = (int)Math.round((value - RESEG_RATIO_BUCKET_MIN) / RESEG_RATIO_BUCKET_STEP);
        index = Math.max(0, Math.min(mBucketCount - 1, index));
        mCounts[index] += weight;
    }

    public List<Bucket> buildSmoothedSeries()
    {
        // rolling-average, edge-trimmed (buckets without a full window are dropped, then duplicates are removed
        int half = RESEG_ROLLING_AVG_WINDOW / 2;
        List<Bucket> smoothed = Lists.newArrayList();

        for(int i = half; i < mBucketCount - half; i++)
        {
            double sum = 0;
            for(int k = i - half; k <= i + half; k++)
            {
                sum += mCounts[k];
            }

            double avg = sum / RESEG_ROLLING_AVG_WINDOW;
            double level = round(RESEG_RATIO_BUCKET_MIN + i * RESEG_RATIO_BUCKET_STEP, 2);
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
