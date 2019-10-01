package com.hartwig.hmftools.sage.context;

@Deprecated
public class ReadContextCounterOld extends ReadContext {

    private int count;

    public ReadContextCounterOld(final ReadContext readContext) {
        super(readContext.readBypePosition, readContext.bytes);
        this.count = 0;
    }

    public void increment() {
        count++;
    }

    public int count() {
        return count;
    }

}
