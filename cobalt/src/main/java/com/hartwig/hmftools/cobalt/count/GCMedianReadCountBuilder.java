package com.hartwig.hmftools.cobalt.count;

import java.util.HashMap;
import java.util.Map;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.genome.gc.GCBucket;
import com.hartwig.hmftools.common.genome.gc.GCMedianReadCount;
import com.hartwig.hmftools.common.genome.gc.GCProfile;

import org.jetbrains.annotations.NotNull;

public class GCMedianReadCountBuilder {

    private static final int MIN_BUCKET = 20;
    private static final int MAX_BUCKET = 60;

    private final boolean useInterpolatedMedian;

    private final ReadCountMedian medianSample;
    private final Map<GCBucket, ReadCountMedian> medianPerGCBucket;

    public GCMedianReadCountBuilder(boolean useInterpolatedMedian) {
        this.useInterpolatedMedian = useInterpolatedMedian;
        medianSample = new ReadCountMedian();
        medianPerGCBucket = new HashMap<>();
    }

    public void add(@NotNull final GCProfile profile, @NotNull final ReadCount readCount) {
        Preconditions.checkArgument(profile.contains(readCount));
        add(profile, readCount.readCount());
    }

    public void add(@NotNull final GCProfile profile, double readCount) {
        final GCBucket gcBucket = GCBucket.create(profile);

        if (gcBucket.bucket() >= MIN_BUCKET && gcBucket.bucket() <= MAX_BUCKET) {
            medianSample.addRead(readCount);
            medianPerGCBucket.computeIfAbsent(gcBucket, integer -> new ReadCountMedian()).addRead(readCount);
        }
    }

    @NotNull
    public GCMedianReadCount build() {
        final Map<GCBucket, Double> gcBucketMeans = new HashMap<>();
        for (GCBucket gcBucket : medianPerGCBucket.keySet()) {
            ReadCountMedian medianCalc = medianPerGCBucket.get(gcBucket);
            gcBucketMeans.put(gcBucket, useInterpolatedMedian ? medianCalc.interpolatedMedian() : medianCalc.median());
        }

        return new GCMedianReadCount(medianSample.mean(),
                useInterpolatedMedian ? medianSample.interpolatedMedian() : medianSample.median(),
                gcBucketMeans);
    }
}

