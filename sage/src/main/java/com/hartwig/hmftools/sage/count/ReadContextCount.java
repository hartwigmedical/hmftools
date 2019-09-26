package com.hartwig.hmftools.sage.count;

public class ReadContextCount extends ReadContext {

    private int count;

    public ReadContextCount(final ReadContext readContext) {
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
