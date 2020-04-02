package com.hartwig.hmftools.sage.quality;

public class QualityCount {
    private final QualityRecord key;
    private int count;

    public QualityCount(final QualityRecord key) {
        this.key = key;
    }

    public int count() {
        return count;
    }

    public void increment() {
        count++;
    }

    public void increment(int increment) {
        count += increment;
    }

    int position() {
        return key.position();
    }

    byte ref() {
        return key.ref();
    }

    byte alt() {
        return key.alt();
    }

    byte qual() {
        return key.qual();
    }

    boolean firstOfPair() {
        return key.firstOfPair();
    }

}
