package com.hartwig.hmftools.sage.count;

import static com.hartwig.hmftools.sage.count.ReadContext.ReadContextMatch.FULL;
import static com.hartwig.hmftools.sage.count.ReadContext.ReadContextMatch.NONE;
import static com.hartwig.hmftools.sage.count.ReadContext.ReadContextMatch.PARTIAL;

import java.util.Arrays;

import org.jetbrains.annotations.NotNull;

public class ReadContext {

    enum ReadContextMatch {
        NONE,
        PARTIAL,
        FULL
    }

    private static final int LENGTH = 10;

    final byte[] bytes;
    final int readBypePosition;

    private final String alt;

    public ReadContext(int readBytePosition, byte[] readBytes) {
        final byte readByte = readBytes[readBytePosition];
        alt = String.valueOf((char) readByte);

        this.bytes = readBytes;
        this.readBypePosition = readBytePosition;
    }

    private int minIndex() {
        return Math.max(0, readBypePosition - LENGTH);
    }

    private int maxIndex() {
        return Math.min(bytes.length - 1, readBypePosition + LENGTH);
    }

    private int rightLength() {
        return maxIndex() - readBypePosition;
    }

    private int leftLength() {
        return readBypePosition - minIndex();
    }

    public String alt() {
        return alt;
    }

    public boolean isComplete() {
        return maxIndex() - minIndex() == 2 * LENGTH;
    }

    public ReadContextMatch match(@NotNull final ReadContext other) {

        if (!isComplete() && !other.isComplete()) {
            return NONE;
        }

        for (int i = 0; i <= Math.min(rightLength(), other.rightLength()); i++) {
            if (bytes[readBypePosition + i] != other.bytes[other.readBypePosition + i]) {
                return NONE;
            }
        }

        for (int i = 1; i <= Math.min(leftLength(), other.leftLength()); i++) {
            if (bytes[readBypePosition - i] != other.bytes[other.readBypePosition - i]) {
                return NONE;
            }
        }

        return isComplete() && other.isComplete() ? FULL : PARTIAL;
    }

    @Override
    public String toString() {
        return new String(Arrays.copyOfRange(bytes, minIndex(), maxIndex() + 1));
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof ReadContext)) {
            return false;
        }
        final ReadContext that = (ReadContext) o;
        return match(that) == FULL;
    }

    @Override
    public int hashCode() {
        int result = 1;

        for (int i = minIndex(); i <= maxIndex(); i++) {
            result = 31 * result + bytes[i];
        }

        return result;
    }
}
