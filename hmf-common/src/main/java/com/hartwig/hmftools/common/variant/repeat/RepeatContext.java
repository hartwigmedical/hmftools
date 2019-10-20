package com.hartwig.hmftools.common.variant.repeat;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class RepeatContext {
    private final byte[] bases;
    private final int position;
    private final int length;
    private final int forwardCount;
    private final int backwardCount;

    public RepeatContext(final byte[] bases, final int startIndex, final int length, int forwardCount, int backwardCount) {
        this.bases = bases;
        this.position = startIndex;
        this.length = length;
        this.backwardCount = backwardCount;
        this.forwardCount = forwardCount;
    }

    public int startIndex() {
        return position - backwardCount * length;
    }

    public int endIndex() {
        return position + forwardCount * length - 1;
    }

    public int count() {
        return forwardCount + backwardCount;
    }

    @NotNull
    public String sequence() {
        return toString();
    }

    @Override
    public String toString() {
        return length > 0 ? new String(bases, position, length) : Strings.EMPTY;
    }
}