package com.hartwig.hmftools.sage.read;

import static com.hartwig.hmftools.sage.read.ReadContextMatch.FULL;
import static com.hartwig.hmftools.sage.read.ReadContextMatch.NONE;
import static com.hartwig.hmftools.sage.read.ReadContextMatch.PARTIAL;

import java.util.Arrays;

import org.jetbrains.annotations.NotNull;

public class IndexedBases {

    @NotNull
    public static IndexedBases resize(final int position, final int recordIndex, final int recordLeftCentreIndex,
            final int recordRightCentreIndex, final int flankSize, final byte[] recordBases) {
        int recordLeftFlankIndex = Math.max(0, recordLeftCentreIndex - flankSize);
        int recordLeftFlankLength = recordLeftCentreIndex - recordLeftFlankIndex;
        int recordRightFlankIndex = Math.min(recordBases.length - 1, recordRightCentreIndex + flankSize);

        int rightCentreIndex = recordLeftFlankLength + recordRightCentreIndex - recordLeftCentreIndex;
        int index = recordLeftFlankLength + recordIndex - recordLeftCentreIndex;
        byte[] bases = Arrays.copyOfRange(recordBases, recordLeftFlankIndex, recordRightFlankIndex + 1);
        return new IndexedBases(position, index, recordLeftFlankLength, rightCentreIndex, flankSize, bases);
    }

    @NotNull
    public static IndexedBases refCentre(@NotNull final IndexedBases readSequences, @NotNull final IndexedBases refSequence) {

        int otherRefIndex = readSequences.position() - refSequence.position() + refSequence.index();
        int otherLeftCentreIndex = readSequences.otherLeftCentreIndex(otherRefIndex);
        int otherRightCentreIndex = readSequences.otherRightCentreIndex(otherRefIndex);

        return resize(readSequences.position(), otherRefIndex, otherLeftCentreIndex, otherRightCentreIndex, 0, refSequence.bases());
    }

    private final int position;
    private final int index;
    private final int flankSize;
    private final int leftFlankIndex;
    private final int leftCentreIndex;
    private final int rightCentreIndex;
    private final int rightFlankIndex;
    private final byte[] bases;

    public IndexedBases(final int position, final int index, final byte[] bases) {
        this.position = position;
        this.index = index;
        this.leftCentreIndex = index;
        this.rightCentreIndex = index;
        this.bases = bases;
        this.leftFlankIndex = index;
        this.rightFlankIndex = index;
        this.flankSize = 0;
    }

    public IndexedBases(final int position, final int index, final int leftCentreIndex, int rightCentreIndex, int flankSize,
            final byte[] bases) {
        this.position = position;
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

        int leftFlankingBases = leftFlankMatchingBases(otherReadIndex, other.bases);
        if (leftFlankingBases < 0) {
            return false;
        }

        int rightFlankingBases = rightFlankMatchingBases(otherReadIndex, other.bases);
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

    @NotNull
    public ReadContextMatch matchAtPosition(int otherReadIndex, byte[] otherBases) {

        if (otherReadIndex < 0 || !flanksComplete()) {
            return NONE;
        }

        boolean centreMatch = centreMatch(otherReadIndex, otherBases);
        if (!centreMatch) {
            return NONE;
        }

        int leftFlankingBases = leftFlankMatchingBases(otherReadIndex, otherBases);
        if (leftFlankingBases < 0) {
            return NONE;
        }

        int rightFlankingBases = rightFlankMatchingBases(otherReadIndex, otherBases);
        if (rightFlankingBases < 0) {
            return NONE;
        }

        if (leftFlankingBases != flankSize && rightFlankingBases != flankSize) {
            return NONE;
        }

        return leftFlankingBases == rightFlankingBases ? FULL : PARTIAL;
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

    int rightFlankMatchingBases(int otherRefIndex, byte[] otherBases) {
        int otherRightCentreIndex = otherRefIndex + rightCentreIndex - index;
        int otherRightFlankLength = Math.min(otherBases.length - 1, otherRightCentreIndex + this.flankSize) - otherRightCentreIndex;
        int maxLength = Math.min(rightFlankLength(), otherRightFlankLength);

        for (int i = 1; i <= maxLength; i++) {
            if (bases[rightCentreIndex + i] != otherBases[otherRightCentreIndex + i]) {
                return -1;
            }
        }

        return maxLength;
    }

    int leftFlankMatchingBases(int otherRefIndex, byte[] otherBases) {
        int otherLeftCentreIndex = otherRefIndex + leftCentreIndex - index;
        int otherLeftFlankLength = otherLeftCentreIndex - Math.max(0, otherLeftCentreIndex - this.flankSize);
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

    public int flankSize() {
        return flankSize;
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

    public int leftCentreIndex() {
        return leftCentreIndex;
    }

    public int rightCentreIndex() {
        return rightCentreIndex;
    }

    public int position() {
        return position;
    }

}
