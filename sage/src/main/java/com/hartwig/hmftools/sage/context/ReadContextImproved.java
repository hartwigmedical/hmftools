package com.hartwig.hmftools.sage.context;

import org.jetbrains.annotations.NotNull;

public class ReadContextImproved {

    private final int position;
    private final int leftIndex;
    private final int rightIndex;
    private final int buffer;
    private final byte[] bases;

    public ReadContextImproved(final int position, final int leftIndex, final int rightIndex, final int buffer, final byte[] bases) {
        assert (leftIndex > 0);
        assert (rightIndex < bases.length);
        assert (rightIndex >= leftIndex);

        this.position = position;
        this.buffer = buffer;
        this.leftIndex = leftIndex;
        this.rightIndex = rightIndex;
        this.bases = bases;
    }

    public int position() {
        return position;
    }

    public boolean isComplete() {
        return leftFlankLength() == buffer && rightFlankLength() == buffer;
    }

    @NotNull
    public ReadContextMatch matchAtPosition(@NotNull final ReadContextImproved other) {

        if (centreLength() != other.centreLength() || position != other.position) {
            return ReadContextMatch.NONE;
        }

        if (!centerMatch(other) || !leftFlankMatch(other) || !rightFlankMatch(other)) {
            return ReadContextMatch.NONE;
        }

        return isComplete() && other.isComplete() ? ReadContextMatch.FULL : ReadContextMatch.PARTIAL;
    }

    public boolean isWithin(byte[] bases) {
        if (!isComplete()) {
            return false;
        }

        for (int i = 0; i < bases.length - length() - 1; i++) {
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

            if (bases[index + i] != this.bases[startIndex + i]) {
                return false;
            }
        }
        return true;
    }

    private boolean centerMatch(@NotNull final ReadContextImproved other) {
        assert other.centreLength() == centreLength();

        for (int i = 0; i < centreLength(); i++) {
            if (bases[leftIndex + i] != other.bases[other.leftIndex + i]) {
                return false;
            }
        }

        return true;
    }

    private boolean leftFlankMatch(@NotNull final ReadContextImproved other) {

        for (int i = 1; i <= Math.min(leftFlankLength(), other.leftFlankLength()); i++) {
            if (bases[leftIndex - i] != other.bases[other.leftIndex - i]) {
                return false;
            }
        }

        return true;
    }

    private boolean rightFlankMatch(@NotNull final ReadContextImproved other) {
        assert other.rightFlankLength() == rightFlankLength();

        for (int i = 1; i <= Math.min(rightFlankLength(), other.rightFlankLength()); i++) {
            if (bases[rightIndex + i] != other.bases[other.rightIndex + i]) {
                return false;
            }
        }

        return true;
    }

    private int leftFlankStartIndex() {
        return Math.max(0, leftIndex - buffer);
    }

    private int leftFlankLength() {
        return leftIndex - leftFlankStartIndex();
    }

    private int rightFlankEndIndex() {
        return Math.min(bases.length - 1, rightIndex + buffer);
    }

    private int rightFlankLength() {
        return rightFlankEndIndex() - rightIndex;
    }

    private int centreLength() {
        return rightIndex - leftIndex + 1;
    }

    private int length() {
        return rightFlankEndIndex() - leftFlankStartIndex() + 1;
    }

    @Override
    public String toString() {
        return new String(bases, leftFlankStartIndex(), length());
    }

}
