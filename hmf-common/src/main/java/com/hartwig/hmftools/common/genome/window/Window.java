package com.hartwig.hmftools.common.genome.window;

public class Window {

    private final int mSize;

    public Window(final int size) {
        this.mSize = size;
    }

    public int getSize() { return mSize; }

    public int start(int position) {
        return (position - 1) / mSize * mSize + 1;
    }

    public int end(int position) {
        return start(position) + mSize - 1;
    }
}
