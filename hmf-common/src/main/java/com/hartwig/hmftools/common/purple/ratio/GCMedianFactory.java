package com.hartwig.hmftools.common.purple.ratio;

import java.util.Map;

import com.google.common.collect.Maps;

class GCMedianFactory {

    private static final int MIN_BUCKET = 20;
    private static final int MAX_BUCKET = 60;

    private final IntegerMedian global;
    private final Map<Integer, IntegerMedian> gcContentMedian;

    GCMedianFactory() {
        gcContentMedian = Maps.newHashMap();
        global = new IntegerMedian();
    }

    void addReadCount(int gcBucket, int readCount) {
        assert (gcBucket <= 100);
        if (gcBucket >= MIN_BUCKET && gcBucket <= MAX_BUCKET) {
            global.addRead(readCount);
            gcContentMedian.computeIfAbsent(gcBucket, integer -> new IntegerMedian()).addRead(readCount);
        }
    }

    double globalNormalFactor() {
        return 1.0 * global.median() / global.mean();
    }

    Map<Integer, GCMedian> medianCountPerGCBucket() {

        final Map<Integer, GCMedian> result = Maps.newHashMap();

        for (Integer key : gcContentMedian.keySet()) {
            result.put(key, ImmutableGCMedian.builder().gcContent(key).medianCount(gcContentMedian.get(key).median()).build());
        }

        return result;
    }

}
