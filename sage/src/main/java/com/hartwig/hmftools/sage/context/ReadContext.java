package com.hartwig.hmftools.sage.context;

import static com.hartwig.hmftools.sage.context.ReadContext.ReadContextMatch.FULL;
import static com.hartwig.hmftools.sage.context.ReadContext.ReadContextMatch.NONE;
import static com.hartwig.hmftools.sage.context.ReadContext.ReadContextMatch.PARTIAL;

import java.util.Arrays;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class ReadContext {

    enum ReadContextMatch {
        NONE,
        PARTIAL,
        FULL
    }

    private static final int DEFAULT_BUFFER = 20;

    private final byte[] readBytes;
    private final int readBytePosition;
    private final int buffer;

    private final String alt;
    private final int distance;
    private final String difference;

    public ReadContext(int readBytePosition, byte[] readBytes) {
        this(DEFAULT_BUFFER, readBytePosition, readBytes);
    }

    public ReadContext(int readBytePosition, @NotNull final SAMRecord record, byte[] refBases) {
        this(DEFAULT_BUFFER, readBytePosition, record, refBases);
    }

    public ReadContext(int buffer, int readBytePosition, @NotNull final SAMRecord record, byte[] refBases) {
        this.buffer = buffer;
        this.readBytes = record.getReadBases();
        this.readBytePosition = readBytePosition;

        final ReadContextDistance readContextDistance = new ReadContextDistance(buffer, readBytePosition, record, refBases);
        this.distance = readContextDistance.distance();
        this.difference = readContextDistance.cigar();
        this.alt = String.valueOf((char) readBytes[readBytePosition]);
    }

    ReadContext(int buffer, int readBytePosition, byte[] readBytes) {
        alt = String.valueOf((char) readBytes[readBytePosition]);

        this.buffer = buffer;
        this.readBytes = readBytes;
        this.readBytePosition = readBytePosition;

        this.difference = Strings.EMPTY;
        this.distance = -1;
    }

    private int minIndex() {
        return Math.max(0, readBytePosition - buffer);
    }

    private int rightLength(final int readBytePosition, final byte[] bytes) {
        return maxIndex(readBytePosition, bytes) - readBytePosition;
    }

    private int leftLength(final int readBytePosition) {
        return readBytePosition - minIndex(readBytePosition);
    }

    public byte[] readBytes() {
        return readBytes;
    }

    public int readBytePosition() {
        return readBytePosition;
    }

    private int maxIndex() {
        return Math.min(readBytes.length - 1, readBytePosition + buffer);
    }

    public String alt() {
        return alt;
    }

    public boolean isComplete() {
        return maxIndex() - minIndex() == 2 * buffer;
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

    private int minIndex(final int readBypePosition) {
        return Math.max(0, readBypePosition - buffer);
    }

    private int maxIndex(final int readBytePosition, final byte[] bytes) {
        return Math.min(bytes.length - 1, readBytePosition + buffer);
    }

    private boolean isComplete(final int readBytePosition, final byte[] bytes) {
        return maxIndex(readBytePosition, bytes) - minIndex(readBytePosition) == 2 * buffer;
    }

    public int distanceFromRef() {
        return distance;
    }

    public String differenceFromRef() {
        return difference;
    }
}
