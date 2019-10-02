package com.hartwig.hmftools.sage.context;

import static com.hartwig.hmftools.sage.context.ReadContext.ReadContextMatch.FULL;
import static com.hartwig.hmftools.sage.context.ReadContext.ReadContextMatch.NONE;
import static com.hartwig.hmftools.sage.context.ReadContext.ReadContextMatch.PARTIAL;

import java.util.Arrays;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class ReadContext {

    enum ReadContextMatch {
        NONE,
        PARTIAL,
        FULL
    }

    private static final int LENGTH = 15;

    private final byte[] readBytes;
    private final int readBytePosition;

    private final String alt;
    private final int distance;
    private final String difference;

    public ReadContext(int readBytePosition, byte[] readBytes, int refBytePosition, byte[] refBytes) {
        alt = String.valueOf((char) readBytes[readBytePosition]);

        this.readBytes = readBytes;
        this.readBytePosition = readBytePosition;

        if (isComplete(readBytePosition, readBytes) && isComplete(refBytePosition, refBytes)) {

            int length = 2 * LENGTH + 1;
            int refStartIndex = minIndex(refBytePosition);
            int readStartIndex = minIndex(readBytePosition);

            int distance = 0;
            char[] diffChar = new char[length];
            for (int i = 0; i < length; i++) {
                byte refByte = refBytes[refStartIndex + i];
                byte readByte = readBytes[readStartIndex + i];

                if (refByte == readByte) {
                    diffChar[i] = '.';
                } else {
                    diffChar[i] = 'X';
                    distance++;
                }
            }

            this.distance = distance;
            difference = new String(diffChar);

        } else {
            this.difference = Strings.EMPTY;
            this.distance = -1;
        }
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
        return readBytes;
    }

    public int readBytePosition() {
        return readBytePosition;
    }

    private int maxIndex() {
        return Math.min(readBytes.length - 1, readBytePosition + LENGTH);
    }

    public String alt() {
        return alt;
    }

    public boolean isComplete() {
        return maxIndex() - minIndex() == 2 * LENGTH;
    }

    public ReadContextMatch match(@NotNull final ReadContext other) {
        return match(other.readBytePosition, other.readBytes);
    }

    public ReadContextMatch match(int otherReadBytePosition, byte[] otherReadBytes) {

        if (!isComplete() && !isComplete(otherReadBytePosition, otherReadBytes)) {
            return NONE;
        }

        for (int i = 0; i <= Math.min(rightLength(readBytePosition, readBytes), rightLength(otherReadBytePosition, otherReadBytes)); i++) {
            if (readBytes[this.readBytePosition + i] != otherReadBytes[otherReadBytePosition + i]) {
                return NONE;
            }
        }

        for (int i = 1; i <= Math.min(leftLength(readBytePosition), leftLength(otherReadBytePosition)); i++) {
            if (readBytes[this.readBytePosition - i] != otherReadBytes[otherReadBytePosition - i]) {
                return NONE;
            }
        }

        return isComplete() && isComplete(otherReadBytePosition, otherReadBytes) ? FULL : PARTIAL;
    }

    @Override
    public String toString() {
        return new String(Arrays.copyOfRange(readBytes, minIndex(), maxIndex() + 1));
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
            result = 31 * result + readBytes[i];
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

    public int distanceFromRef() {
        return distance;
    }

    public String differenceFromRef() {
        return difference;
    }
}
