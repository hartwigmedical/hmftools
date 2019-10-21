package com.hartwig.hmftools.common.variant;

import org.apache.logging.log4j.util.Strings;

public class MicrohomologyContext {

    private final byte[] bases;
    private final int leftAlignedPosition;
    private final int length;

    public MicrohomologyContext(final int leftAlignedPosition, final byte[] bases, final int length) {
        this.bases = bases;
        this.leftAlignedPosition = leftAlignedPosition;
        this.length = length;
    }

    public int position() {
        return leftAlignedPosition;
    }

    public byte[] readSequence() {
        return bases;
    }

    @Override
    public String toString() {
        return length > 0 ? new String(bases, homologyIndex(), length) : Strings.EMPTY;
    }

    public int length() {
        return length;
    }

    public int homologyIndex() {
        return leftAlignedPosition + 1;
    }
}
