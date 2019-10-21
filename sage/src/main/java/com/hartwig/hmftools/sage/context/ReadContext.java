package com.hartwig.hmftools.sage.context;

import static com.hartwig.hmftools.sage.context.ReadContextMatch.FULL;
import static com.hartwig.hmftools.sage.context.ReadContextMatch.NONE;
import static com.hartwig.hmftools.sage.context.ReadContextMatch.PARTIAL;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class ReadContext {

    private final int position;
    private final int readIndex;
    private final int leftCentreIndex;
    private final int rightCentreIndex;
    private final int flankSize;
    private final byte[] readBases;
    private final int distance;
    private final String distanceCigar;
    private final String repeat;
    private final String microhomology;
    private final int jitter;

    public ReadContext(final String microhomology, final String repeat, final int refPosition, final int readIndex,
            final int leftCentreIndex, final int rightCentreIndex, final int flankSize, final byte[] readBases) {
        assert (leftCentreIndex > 0);
        assert (rightCentreIndex < readBases.length);
        assert (rightCentreIndex >= leftCentreIndex);

        this.position = refPosition;
        this.flankSize = flankSize;
        this.leftCentreIndex = leftCentreIndex;
        this.rightCentreIndex = rightCentreIndex;
        this.readBases = readBases;
        this.readIndex = readIndex;
        this.distance = 0;
        this.distanceCigar = Strings.EMPTY;
        this.jitter = repeat.length() >= centreLength() ? 0 : repeat.length();
        this.repeat = repeat;
        this.microhomology = microhomology;
    }

    public ReadContext(final String microhomology, final String repeat, final int refPosition, final int readIndex,
            final int leftCentreIndex, final int rightCentreIndex, final int flankSize, byte[] refBases, @NotNull final SAMRecord record) {
        assert (leftCentreIndex >= 0);
        assert (rightCentreIndex >= leftCentreIndex);

        this.position = refPosition;
        this.flankSize = flankSize;
        this.leftCentreIndex = leftCentreIndex;
        this.rightCentreIndex = rightCentreIndex;
        this.readBases = record.getReadBases();
        this.readIndex = readIndex;

        ReadContextDistance distance = new ReadContextDistance(leftFlankStartIndex(), rightFlankEndIndex(), record, refBases);
        this.distance = distance.distance();
        this.distanceCigar = distance.cigar();
        this.jitter = repeat.length() >= centreLength() ? 0 : repeat.length();
        this.repeat = repeat;
        this.microhomology = microhomology;
    }

    public int position() {
        return position;
    }

    int distanceFromReadEdge(int readIndex, SAMRecord record) {
        int leftOffset = this.readIndex - leftCentreIndex;
        int rightOffset = rightCentreIndex - this.readIndex;

        int leftIndex = readIndex - leftOffset;
        int rightIndex = readIndex + rightOffset;

        return Math.min(leftIndex, record.getReadBases().length - rightIndex - 1);
    }

    public boolean isComplete() {
        return leftFlankLength() == flankSize && rightFlankLength() == flankSize;
    }

    public boolean isFullMatch(@NotNull final ReadContext other) {
        return isComplete() && other.isComplete() && centreMatch(other.readIndex, other.readBases) == ReadContextMatch.FULL
                && leftFlankMatchingBases(other.readIndex, other.readBases) == flankSize
                && rightFlankMatchingBases(other.readIndex, other.readBases) == flankSize;
    }

    public boolean phased(@NotNull final ReadContext other) {
        int offset = position - other.position;
        int otherReadIndex = other.readIndex + offset;

        ReadContextMatch centreMatch = centreMatch(otherReadIndex, other.readBases);
        if (centreMatch != FULL) {
            return false;
        }

        final int rightFlankSize;
        final int leftFlankSize;
        if (offset < 0) {
            leftFlankSize = flankSize + offset;
            rightFlankSize = flankSize;
        } else {
            leftFlankSize = flankSize;
            rightFlankSize = flankSize - offset;
        }

        int leftFlankingBases = leftFlankMatchingBases(otherReadIndex, other.readBases, leftFlankSize);
        if (leftFlankingBases < 0) {
            return false;
        }

        int rightFlankingBases = rightFlankMatchingBases(otherReadIndex, other.readBases, rightFlankSize);
        return rightFlankingBases >= 0;
    }

    @NotNull
    public ReadContextMatch matchAtPosition(int otherReadIndex, byte[] otherBases) {

        if (!isComplete()) {
            return NONE;
        }

        ReadContextMatch centreMatch = centreMatch(otherReadIndex, otherBases);
        if (centreMatch == NONE) {
            return NONE;
        }

        int leftFlankingBases = leftFlankMatchingBases(otherReadIndex, otherBases);
        if (leftFlankingBases < 0) {
            return NONE;
        }

        final int rightFlankingStartIndex;
        switch (centreMatch) {
            case FULL:
                rightFlankingStartIndex = otherReadIndex;
                break;
            case JITTER_ADDED:
                rightFlankingStartIndex = otherReadIndex + jitter;
                break;
            case JITTER_REMOVED:
                rightFlankingStartIndex = otherReadIndex - jitter;
                break;
            default:
                throw new IllegalStateException("Unable to handle centre match type " + centreMatch);
        }

        int rightFlankingBases = rightFlankMatchingBases(rightFlankingStartIndex, otherBases);
        if (rightFlankingBases < 0) {
            return NONE;
        }

        if (leftFlankingBases != flankSize && rightFlankingBases != flankSize) {
            return NONE;
        }

        if (centreMatch == ReadContextMatch.FULL) {
            return leftFlankingBases == rightFlankingBases ? FULL : PARTIAL;
        }

        return centreMatch;
    }

    public boolean isWithin(byte[] bases) {
        if (!isComplete()) {
            return false;
        }

        for (int i = 0; i < bases.length - length(); i++) {
            if (isWithin(i, bases)) {
                return true;
            }
        }

        return false;
    }

    private boolean isWithin(int index, byte[] bases) {
        int startIndex = leftFlankStartIndex();
        for (int i = 0; i < length(); i++) {
            if (index + i >= bases.length) {
                return false;
            }

            if (bases[index + i] != this.readBases[startIndex + i]) {
                return false;
            }
        }
        return true;
    }

    @VisibleForTesting
    @NotNull
    ReadContextMatch centreMatch(int otherRefIndex, byte[] otherBases) {
        int otherLeftCentreIndex = otherRefIndex + leftCentreIndex - readIndex;
        if (otherLeftCentreIndex < 0) {
            return NONE;
        }

        // Up to jitter
        int lengthUntilJitter = centreLength() - 1 - jitter;
        if (otherLeftCentreIndex + lengthUntilJitter - 1 >= otherBases.length) {
            return NONE;
        }
        for (int i = 0; i < lengthUntilJitter; i++) {
            if (readBases[leftCentreIndex + i] != otherBases[otherLeftCentreIndex + i]) {
                return NONE;
            }
        }

        boolean fullMatch = true;
        for (int i = lengthUntilJitter; i < centreLength(); i++) {
            if (otherLeftCentreIndex + i >= otherBases.length || readBases[leftCentreIndex + i] != otherBases[otherLeftCentreIndex + i]) {
                fullMatch = false;
                break;
            }
        }

        if (fullMatch) {
            return FULL;
        }

        if (jitter == 0) {
            return NONE;
        }

        if (otherLeftCentreIndex + lengthUntilJitter < otherBases.length && readBases[rightCentreIndex] == otherBases[otherLeftCentreIndex
                + lengthUntilJitter]) {
            return ReadContextMatch.JITTER_REMOVED;
        }

        int otherRightCentreIndex = otherLeftCentreIndex + centreLength() - 1 + jitter;
        for (int i = 0; i <= jitter; i++) {
            if (otherRightCentreIndex - jitter + i >= otherBases.length || readBases[rightCentreIndex - jitter + i] != otherBases[
                    otherRightCentreIndex - jitter + i]) {
                return ReadContextMatch.NONE;
            }
        }

        return ReadContextMatch.JITTER_ADDED;
    }

    @VisibleForTesting
    int leftFlankMatchingBases(int otherRefIndex, byte[] otherBases) {
        return leftFlankMatchingBases(otherRefIndex, otherBases, flankSize);
    }

    int leftFlankMatchingBases(int otherRefIndex, byte[] otherBases, int flankSize) {
        int otherLeftCentreIndex = otherRefIndex + leftCentreIndex - readIndex;
        int otherLeftFlankLength = otherLeftCentreIndex - Math.max(0, otherLeftCentreIndex - flankSize);
        int totalLength = Math.min(leftFlankLength(), otherLeftFlankLength);

        for (int i = 1; i <= totalLength; i++) {
            if (readBases[leftCentreIndex - i] != otherBases[otherLeftCentreIndex - i]) {
                return -1;
            }
        }

        return totalLength;
    }

    @VisibleForTesting
    int rightFlankMatchingBases(int otherRefIndex, byte[] otherBases) {
        return rightFlankMatchingBases(otherRefIndex, otherBases, flankSize);
    }

    int rightFlankMatchingBases(int otherRefIndex, byte[] otherBases, int flankSize) {
        int otherRightCentreIndex = otherRefIndex + rightCentreIndex - readIndex;
        int otherRightFlankLength = Math.min(otherBases.length - 1, otherRightCentreIndex + flankSize) - otherRightCentreIndex;
        int maxLength = Math.min(rightFlankLength(), otherRightFlankLength);

        for (int i = 1; i <= maxLength; i++) {
            if (readBases[rightCentreIndex + i] != otherBases[otherRightCentreIndex + i]) {
                return -1;
            }
        }

        return maxLength;
    }

    int centreStartIndex() {
        return leftCentreIndex;
    }

    int centreEndIndex() {
        return rightCentreIndex;
    }

    private int leftFlankStartIndex() {
        return Math.max(0, leftCentreIndex - flankSize);
    }

    private int leftFlankLength() {
        return leftCentreIndex - leftFlankStartIndex();
    }

    private int rightFlankEndIndex() {
        return Math.min(readBases.length - 1, rightCentreIndex + flankSize);
    }

    private int rightFlankLength() {
        return rightFlankEndIndex() - rightCentreIndex;
    }

    private int centreLength() {
        return rightCentreIndex - leftCentreIndex + 1;
    }

    private int length() {
        return rightFlankEndIndex() - leftFlankStartIndex() + 1;
    }

    @Override
    public String toString() {
        return new String(readBases, leftFlankStartIndex(), length());
    }

    public int distance() {
        return distance;
    }

    public String distanceCigar() {
        return distanceCigar;
    }

    @VisibleForTesting
    @NotNull
    String centerBases() {
        return new String(readBases, leftCentreIndex, rightCentreIndex - leftCentreIndex + 1);
    }

    @NotNull
    public String microhomology() {
        return microhomology;
    }

    @NotNull
    public String repeat() {
        return repeat;
    }
}
