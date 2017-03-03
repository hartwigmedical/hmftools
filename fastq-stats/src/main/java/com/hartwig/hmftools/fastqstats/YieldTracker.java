package com.hartwig.hmftools.fastqstats;

public class YieldTracker implements Tracker{
    private int count = 0;

    @Override
    public void addValue(final int v) {
        count ++;
    }

    @Override
    public long getCount() {
        return count;
    }
}
