package com.hartwig.hmftools.window;

public class Window {

    private final int size;

    public Window(final int size) {
        this.size = size;
    }

    public long start(long position) {
        return (position - 1) / size * size + 1;
    }

    public long end(long position) {
        return start(position) + size - 1;
    }
}
