package com.hartwig.hmftools.common.gc;

import java.util.Map;

import org.jetbrains.annotations.NotNull;

class GCMedianReadCountImpl implements GCMedianReadCount {

    private final int mean;
    private final int median;
    private final Map<GCBucket, Integer> medianReadCountPerGCBucket;

    GCMedianReadCountImpl(final int mean, final int median, final Map<GCBucket, Integer> medianReadCountPerGCBucket) {
        this.mean = mean;
        this.median = median;
        this.medianReadCountPerGCBucket = medianReadCountPerGCBucket;
    }

    @Override
    public int meanReadCount() {
        return mean;
    }

    @Override
    public int medianReadCount() {
        return median;
    }

    @Override
    public int medianReadCount(@NotNull final GCBucket bucket) {
        return medianReadCountPerGCBucket.getOrDefault(bucket, -1);
    }

}
