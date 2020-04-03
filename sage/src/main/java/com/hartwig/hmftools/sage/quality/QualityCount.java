package com.hartwig.hmftools.sage.quality;

public class QualityCount implements QualityRecord {
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

    public int position() {
        return key.position();
    }

    public byte ref() {
        return key.ref();
    }

    public byte alt() {
        return key.alt();
    }

    public byte qual() {
        return key.qual();
    }

    @Override
    public int readIndex() {
        return key.readIndex();
    }

    @Override
    public byte[] trinucleotideContext() {
        return key.trinucleotideContext();
    }

    public boolean firstOfPair() {
        return key.firstOfPair();
    }

}
