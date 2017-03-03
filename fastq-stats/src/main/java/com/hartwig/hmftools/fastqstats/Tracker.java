package com.hartwig.hmftools.fastqstats;

public interface Tracker {
    void addValue(int v);

    long getCount();
}
