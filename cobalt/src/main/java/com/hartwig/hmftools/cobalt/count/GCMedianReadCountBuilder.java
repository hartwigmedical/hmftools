package com.hartwig.hmftools.cobalt.count;

import java.util.HashMap;
import java.util.Map;

import com.google.common.base.Preconditions;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.gc.GCBucket;
import com.hartwig.hmftools.common.genome.gc.GCMedianReadCount;
import com.hartwig.hmftools.common.genome.gc.GCProfile;

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

    public void add(@NotNull final GCProfile profile, int readCount) {
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
            gcBucketMeans.put(gcBucket, medianPerGCBucket.get(gcBucket).interpolatedMedian());
        }

        return new GCMedianReadCount(medianSample.mean(), medianSample.interpolatedMedian(), gcBucketMeans);
    }
}

