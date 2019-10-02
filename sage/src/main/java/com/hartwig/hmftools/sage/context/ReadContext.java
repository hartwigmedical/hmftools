package com.hartwig.hmftools.sage.context;

import static com.hartwig.hmftools.sage.context.ReadContext.ReadContextMatch.FULL;
import static com.hartwig.hmftools.sage.context.ReadContext.ReadContextMatch.NONE;
import static com.hartwig.hmftools.sage.context.ReadContext.ReadContextMatch.PARTIAL;

import java.util.Arrays;

import org.jetbrains.annotations.NotNull;

public class ReadContext {

    enum ReadContextMatch {
        NONE,
        PARTIAL,
        FULL
    }

    private static final int LENGTH = 15;

    private final byte[] bytes;
    private final int readBytePosition;

    private final String alt;

    public ReadContext(int readBytePosition, byte[] readBytes) {
        final byte readByte = readBytes[readBytePosition];
        alt = String.valueOf((char) readByte);

        this.bytes = readBytes;
        this.readBytePosition = readBytePosition;
    }

    private int minIndex() {
        return Math.max(0, readBytePosition - LENGTH);
    }

    private static int rightLength(final int readBytePosition, final byte[] bytes) {
        return maxIndex(readBytePosition, bytes) - readBytePosition;
    }

    private static int leftLength(final int readBytePosition) {
        return readBytePosition - minIndex(readBytePosition);
    }

    public byte[] readBytes() {
        return bytes;
    }

    public int readBytePosition() {
        return readBytePosition;
    }

    private int maxIndex() {
        return Math.min(bytes.length - 1, readBytePosition + LENGTH);
    }

    private int rightLength() {
        return maxIndex() - readBytePosition;
    }

    private int leftLength() {
        return readBytePosition - minIndex();
    }

    public String alt() {
        return alt;
    }

    public boolean isComplete() {
        return maxIndex() - minIndex() == 2 * LENGTH;
    }

    public ReadContextMatch match(@NotNull final ReadContext other) {
        return match(other.readBytePosition, other.bytes);
    }

    public ReadContextMatch match(int otherReadBytePosition, byte[] otherReadBytes) {

        if (!isComplete() && !isComplete(otherReadBytePosition, otherReadBytes)) {
            return NONE;
        }

        for (int i = 0; i <= Math.min(rightLength(), rightLength(otherReadBytePosition, otherReadBytes)); i++) {
            if (bytes[this.readBytePosition + i] != otherReadBytes[otherReadBytePosition + i]) {
                return NONE;
            }
        }

        for (int i = 1; i <= Math.min(leftLength(), leftLength(otherReadBytePosition)); i++) {
            if (bytes[this.readBytePosition - i] != otherReadBytes[otherReadBytePosition - i]) {
                return NONE;
            }
        }

        return isComplete() && isComplete(otherReadBytePosition, otherReadBytes) ? FULL : PARTIAL;
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

    private static int minIndex(final int readBypePosition) {
        return Math.max(0, readBypePosition - LENGTH);
    }

    private static int maxIndex(final int readBytePosition, final byte[] bytes) {
        return Math.min(bytes.length - 1, readBytePosition + LENGTH);
    }

    private static boolean isComplete(final int readBytePosition, final byte[] bytes) {
        return maxIndex(readBytePosition, bytes) - minIndex(readBytePosition) == 2 * LENGTH;
    }
}
