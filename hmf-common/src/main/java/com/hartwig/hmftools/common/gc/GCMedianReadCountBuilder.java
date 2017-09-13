package com.hartwig.hmftools.common.gc;

import java.util.Map;

import com.google.common.base.Preconditions;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.cobalt.ReadCount;

import org.jetbrains.annotations.NotNull;

public class GCMedianReadCountBuilder {

    private static final int MIN_BUCKET = 20;
    private static final int MAX_BUCKET = 60;

    private final ReadCountMedian medianSample;
    private final Map<GCBucket, ReadCountMedian> medianPerGCBucket;

    public GCMedianReadCountBuilder() {
        medianSample = new ReadCountMedian();
        medianPerGCBucket = Maps.newHashMap();
    }

    public void add(@NotNull final GCProfile profile, @NotNull final ReadCount readCount) {
        Preconditions.checkArgument(profile.contains(readCount));
        add(profile, readCount.readCount());
    }

    private void add(@NotNull final GCProfile profile, int readCount) {
        final GCBucket gcBucket = GCBucket.create(profile);

        if (gcBucket.bucket() >= MIN_BUCKET && gcBucket.bucket() <= MAX_BUCKET) {
            medianSample.addRead(readCount);
            medianPerGCBucket.computeIfAbsent(gcBucket, integer -> new ReadCountMedian()).addRead(readCount);
        }
    }

    public GCMedianReadCountImpl build() {

        final Map<GCBucket, Integer> gcBucketMeans = Maps.newHashMap();
        for (GCBucket gcBucket : medianPerGCBucket.keySet()) {
            gcBucketMeans.put(gcBucket, medianPerGCBucket.get(gcBucket).median());
        }

        return new GCMedianReadCountImpl(medianSample.mean(), medianSample.median(), gcBucketMeans);
    }
}

