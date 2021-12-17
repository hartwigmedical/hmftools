package com.hartwig.hmftools.common.genome.window;

public class Window {

    private final int size;

    public Window(final int size) {
        this.size = size;
    }

    public int start(int position) {
        return (position - 1) / size * size + 1;
    }

    public int end(int position) {
        return start(position) + size - 1;
    }
}
