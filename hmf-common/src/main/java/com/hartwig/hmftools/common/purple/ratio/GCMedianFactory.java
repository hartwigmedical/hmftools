package com.hartwig.hmftools.common.purple.ratio;

import java.util.Map;

import com.google.common.collect.Maps;

class GCMedianFactory {

    private static final int MIN_BUCKET = 25;
    private static final int MAX_BUCKET = 75;

    private final Map<Integer, IntegerMedian> gcContentMedian;

    GCMedianFactory() {
        gcContentMedian = Maps.newHashMap();
    }

    void addReadCount(int gcBucket, int readCount) {
        assert (gcBucket <= 100);
        if (gcBucket > 0) {
            gcContentMedian.computeIfAbsent(gcBucket, integer -> new IntegerMedian()).addRead(readCount);
        }
    }

    Map<Integer, GCMedians> medianCountPerGCBucket() {

        final Map<Integer, GCMedians> result = Maps.newHashMap();

        int closestBucketToMin = closestBucketToMin();
        int closestBucketToMax = closestBucketToMax();

        for (Integer key : gcContentMedian.keySet()) {
            int bucket = Math.min(closestBucketToMax, Math.max(closestBucketToMin, key));
            result.put(key, ImmutableGCMedians.builder().gcContent(key).medianCount(gcContentMedian.get(bucket).median()).build());
        }

        return result;
    }

    private int closestBucketToMin() {
        int min = MAX_BUCKET;

        for (Integer bucket : gcContentMedian.keySet()) {
            if (bucket >= MIN_BUCKET && bucket < min) {
                min = bucket;
            }
        }

        return min;
    }

    private int closestBucketToMax() {
        int max = MIN_BUCKET;

        for (Integer bucket : gcContentMedian.keySet()) {
            if (bucket <= MAX_BUCKET && bucket > max) {
                max = bucket;
            }
        }

        return max;
    }
}
