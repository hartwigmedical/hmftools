package com.hartwig.hmftools.fastqstats;

import java.util.function.IntPredicate;

public class PredicateTracker implements Tracker{
    private IntPredicate predicate;
    private long count = 0;

    public PredicateTracker(IntPredicate p) {
        predicate = p;
    }

    @Override
    public void addValue(int v) {
        if(predicate.test(v)){
            count ++;
        }
    }

    @Override
    public long getCount() {
        return count;
    }
}
