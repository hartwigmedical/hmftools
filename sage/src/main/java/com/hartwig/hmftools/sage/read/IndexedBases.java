package com.hartwig.hmftools.sage.read;

import org.jetbrains.annotations.NotNull;

public class IndexedBases {

    private final int index;
    private final int flankSize;
    private final int leftFlankIndex;
    private final int leftCentreIndex;
    private final int rightCentreIndex;
    private final int rightFlankIndex;
    private final byte[] bases;

    public IndexedBases(final int index, final int leftCentreIndex, int rightCentreIndex, final byte[] bases) {
        this.index = index;
        this.leftCentreIndex = leftCentreIndex;
        this.rightCentreIndex = rightCentreIndex;
        this.bases = bases;
        this.leftFlankIndex = leftCentreIndex;
        this.rightFlankIndex = rightCentreIndex;
        this.flankSize = 0;
    }

    public IndexedBases(final int index, final int leftCentreIndex, int rightCentreIndex, int flankSize, final byte[] bases) {
        this.index = index;
        this.leftCentreIndex = leftCentreIndex;
        this.rightCentreIndex = rightCentreIndex;
        this.bases = bases;
        this.leftFlankIndex = Math.max(0, leftCentreIndex - flankSize);
        this.rightFlankIndex = Math.min(bases.length - 1, rightCentreIndex + flankSize);
        this.flankSize = flankSize;
    }

    boolean flanksComplete() {
        return leftFlankLength() == flankSize && rightFlankLength() == flankSize;
    }

    public boolean phased(int offset, @NotNull final IndexedBases other) {
        int otherReadIndex = other.index + offset;

        boolean centreMatch = centreMatch(otherReadIndex, other.bases);
        if (!centreMatch) {
            return false;
        }

        boolean otherCentreMatch = other.centreMatch(index - offset, bases);
        if (!otherCentreMatch) {
            return false;
        }

        int leftFlankingBases = leftFlankMatchingBases(otherReadIndex, other.bases, flankSize);
        if (leftFlankingBases < 0) {
            return false;
        }

        int rightFlankingBases = rightFlankMatchingBases(otherReadIndex, other.bases, flankSize);
        return rightFlankingBases >= 0 && (rightFlankingBases >= flankSize || leftFlankingBases >= flankSize);
    }

    public boolean isCentreCovered(int otherReadIndex, byte[] otherBases) {

        int otherLeftCentreIndex = otherLeftCentreIndex(otherReadIndex);
        if (otherLeftCentreIndex < 0) {
            return false;
        }

        int otherRightCentreIndex = otherRightCentreIndex(otherReadIndex);
        return otherRightCentreIndex < otherBases.length;
    }


    boolean centreMatch(final int otherRefIndex, final byte[] otherBases) {

        int otherLeftCentreIndex = otherLeftCentreIndex(otherRefIndex);
        if (otherLeftCentreIndex < 0) {
            return false;
        }

        int otherRightCentreIndex = otherRightCentreIndex(otherRefIndex);
        if (otherRightCentreIndex >= otherBases.length) {
            return false;
        }

        for (int i = 0; i < centreLength(); i++) {
            if (bases[leftCentreIndex + i] != otherBases[otherLeftCentreIndex + i]) {
                return false;
            }
        }

        return true;
    }

    int rightFlankMatchingBases(int otherRefIndex, byte[] otherBases, int flankSize) {
        int otherRightCentreIndex = otherRefIndex + rightCentreIndex - index;
        int otherRightFlankLength = Math.min(otherBases.length - 1, otherRightCentreIndex + flankSize) - otherRightCentreIndex;
        int maxLength = Math.min(rightFlankLength(), otherRightFlankLength);

        for (int i = 1; i <= maxLength; i++) {
            if (bases[rightCentreIndex + i] != otherBases[otherRightCentreIndex + i]) {
                return -1;
            }
        }

        return maxLength;
    }

    int leftFlankMatchingBases(int otherRefIndex, byte[] otherBases, int flankSize) {
        int otherLeftCentreIndex = otherRefIndex + leftCentreIndex - index;
        int otherLeftFlankLength = otherLeftCentreIndex - Math.max(0, otherLeftCentreIndex - flankSize);
        int totalLength = Math.min(leftFlankLength(), otherLeftFlankLength);

        for (int i = 1; i <= totalLength; i++) {
            if (bases[leftCentreIndex - i] != otherBases[otherLeftCentreIndex - i]) {
                return -1;
            }
        }

        return totalLength;
    }

    private int otherLeftCentreIndex(int otherRefIndex) {
        return otherRefIndex + leftCentreIndex - index;
    }

    private int otherRightCentreIndex(int otherRefIndex) {
        return otherRefIndex + rightCentreIndex - index;
    }

    @NotNull
    public String centerString() {
        return new String(bases, leftCentreIndex, centreLength());
    }

    @Override
    public String toString() {
        return new String(bases, leftFlankIndex, length());
    }

    public int index() {
        return index;
    }

    public byte[] bases() {
        return bases;
    }

    private int length() {
        return rightFlankIndex - leftFlankIndex + 1;
    }

    private int leftFlankLength() {
        return leftCentreIndex - leftFlankIndex;
    }

    private int rightFlankLength() {
        return rightFlankIndex - rightCentreIndex;
    }

    private int centreLength() {
        return rightCentreIndex - leftCentreIndex + 1;
    }
}
