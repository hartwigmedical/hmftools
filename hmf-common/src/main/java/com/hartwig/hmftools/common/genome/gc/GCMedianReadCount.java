package com.hartwig.hmftools.common.genome.gc;

import java.util.Map;

import org.jetbrains.annotations.NotNull;

public class GCMedianReadCount {

    private final double mean;
    private final double median;
    private final Map<GCBucket, Double> medianReadCountPerGCBucket;

    public GCMedianReadCount(final double mean, final double median, final Map<GCBucket, Double> medianReadCountPerGCBucket) {
        this.mean = mean;
        this.median = median;
        this.medianReadCountPerGCBucket = medianReadCountPerGCBucket;
    }

    public double meanReadCount() {
        return mean;
    }

    public double medianReadCount() {
        return median;
    }

    public double medianReadCount(@NotNull final GCBucket bucket) { return medianReadCountPerGCBucket.getOrDefault(bucket, -1.0); }

    public double medianReadCount(@NotNull GCProfile profile) {
        return medianReadCount(GCBucket.create(profile));
    }
}
