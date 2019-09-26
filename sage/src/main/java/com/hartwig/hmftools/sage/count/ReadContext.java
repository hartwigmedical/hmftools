package com.hartwig.hmftools.sage.count;

import java.util.Arrays;
import java.util.Objects;

import org.jetbrains.annotations.NotNull;

public class ReadContext {

    private static final int LENGTH = 10;

    private final String left;
    private final String alt;
    private final String right;

    public ReadContext(int readBytePosition, byte[] readBytes) {
        final byte readByte = readBytes[readBytePosition];
        alt = String.valueOf((char) readByte);

        left = new String(Arrays.copyOfRange(readBytes, Math.max(0, readBytePosition - LENGTH), readBytePosition));
        right = new String(Arrays.copyOfRange(readBytes, readBytePosition + 1, Math.min(readBytes.length, readBytePosition + 1 + LENGTH)));
    }

    ReadContext(final String left, final String alt, final String right) {
        this.left = left;
        this.alt = alt;
        this.right = right;
    }

    public boolean isComplete() {
        return left.length() + alt.length() + right.length() == 2 * LENGTH + 1;
    }

    public boolean match(@NotNull final ReadContext other) {
        if (!alt.equals(other.alt)) {
            return false;
        }

        if (left.length() == other.left.length() && !left.equals(other.left)) {
            return false;
        }

        if (right.length() == other.right.length() && !right.equals(other.right)) {
            return false;
        }

        return left.length() == other.left.length() || right.length() == other.right.length();
    }

    @Override
    public String toString() {
        return left + alt + right;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        final ReadContext that = (ReadContext) o;
        return Objects.equals(left, that.left) && Objects.equals(alt, that.alt) && Objects.equals(right, that.right);
    }

    @Override
    public int hashCode() {

        return Objects.hash(left, alt, right);
    }
}
