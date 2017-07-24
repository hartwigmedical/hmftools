package com.hartwig.hmftools.common.purple.ratio;

import java.util.Map;

import com.google.common.collect.Maps;

public class GCMedian {

    private static final int MIN_GC_BUCKET = 25;
    private static final int MAX_GC_BUCKET = 75;

    private final IntegerMedian globalMedian;
    private final Map<Integer, IntegerMedian> gcContentMedian;
    private final Map<Integer, Integer> gcContentMedianValues;

    public GCMedian() {
        globalMedian = new IntegerMedian();
        gcContentMedian = Maps.newHashMap();
        gcContentMedianValues = Maps.newHashMap();
    }

    public void addRead(int gcContent, int read) {
        assert (gcContent >= 0);
        assert (gcContent <= 100);

        globalMedian.addRead(read);
        gcContentMedian.computeIfAbsent(gcContent, integer -> new IntegerMedian()).addRead(read);
    }

    public int globalMedian() {
        return globalMedian.median();
    }

    public int gcContentMedian(int gcContent) {
        return gcContentMedianValues.computeIfAbsent(gcContent, key -> gcContentMedian.get(key).median());
    }

    public Map<Integer, Integer> medianCountPerGCBucket() {

        final Map<Integer, Integer> result = Maps.newHashMap();

        int closestBucketToMin = closestBucketToMin();
        int closestBucketToMax = closestBucketToMax();

        for (Integer key : gcContentMedian.keySet()) {
            int gcBucket = Math.min(closestBucketToMax, Math.max(closestBucketToMin, key));
            result.put(key, gcContentMedian.get(gcBucket).median());
        }

        return result;
    }

    private int closestBucketToMin() {
        int min = 100;

        for (Integer gcBucket : gcContentMedian.keySet()) {
            if (gcBucket >= MIN_GC_BUCKET && gcBucket < min) {
                min = gcBucket;
            }
        }

        return min;
    }

    private int closestBucketToMax() {
        int max = 0;

        for (Integer gcBucket : gcContentMedian.keySet()) {
            if (gcBucket <= MAX_GC_BUCKET && gcBucket > max) {
                max = gcBucket;
            }
        }

        return max;
    }

}
