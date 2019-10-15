package com.hartwig.hmftools.sage.context;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class ReadContextImproved {

    private final int position;
    private final int readIndex;
    private final int leftCentreIndex;
    private final int rightCentreIndex;
    private final int flankSize;
    private final byte[] readBases;
    private final int distance;
    private final String distanceCigar;

    public ReadContextImproved(final int refPosition, final int readIndex, final int leftCentreIndex, final int rightCentreIndex,
            final int flankSize, final byte[] readBases) {
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
    }

    public ReadContextImproved(final int refPosition, final int readIndex, final int leftCentreIndex, final int rightCentreIndex,
            final int flankSize, byte[] refBases, @NotNull final SAMRecord record) {
        assert (leftCentreIndex > 0);
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

    }

    public int position() {
        return position;
    }

    int baseQuality(int readIndex, SAMRecord record) {
        int leftOffset = this.readIndex - leftCentreIndex;
        int rightOffset = rightCentreIndex - this.readIndex;

        int leftIndex = readIndex - leftOffset;
        int rightIndex = readIndex + rightOffset;

        int leftQuality = record.getBaseQualities()[leftIndex];
        if (leftIndex == rightIndex) {
            return leftQuality;
        }

        int rightQuality = record.getBaseQualities()[rightIndex];
        return Math.min(leftQuality, rightQuality);
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

    public boolean isFullMatch(@NotNull final ReadContextImproved other) {
        return isComplete() && other.isComplete() && centerMatch(other.readIndex, other.readBases)
                && leftFlankMatchingBases(other.readIndex, other.readBases) == flankSize
                && rightFlankMatchingBases(other.readIndex, other.readBases) == flankSize;
    }

    @NotNull
    public ReadContextMatch matchAtPosition(int otherReadIndex, byte[] otherBases) {

        if (!isComplete() || !centerMatch(otherReadIndex, otherBases)) {
            return ReadContextMatch.NONE;
        }

        int leftFlankingBases = leftFlankMatchingBases(otherReadIndex, otherBases);
        if (leftFlankingBases < 0) {
            return ReadContextMatch.NONE;
        }

        int rightFlankingBases = rightFlankMatchingBases(otherReadIndex, otherBases);
        if (rightFlankingBases < 0) {
            return ReadContextMatch.NONE;
        }

        return leftFlankingBases == flankSize && rightFlankingBases == flankSize ? ReadContextMatch.FULL : ReadContextMatch.PARTIAL;
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
    boolean centerMatch(int otherRefIndex, byte[] otherBases) {
        int otherLeftCentreIndex = otherRefIndex + leftCentreIndex - readIndex;
        if (otherLeftCentreIndex < 0) {
            return false;
        }

        int otherRightCentreIndex = otherLeftCentreIndex + centreLength() - 1;
        if (otherRightCentreIndex >= otherBases.length) {
            return false;
        }

        for (int i = 0; i < centreLength(); i++) {
            if (readBases[leftCentreIndex + i] != otherBases[otherLeftCentreIndex + i]) {
                return false;
            }
        }

        return true;
    }

    @VisibleForTesting
    int leftFlankMatchingBases(int otherRefIndex, byte[] otherBases) {
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
}
