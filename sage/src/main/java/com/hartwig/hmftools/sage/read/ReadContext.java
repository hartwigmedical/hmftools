package com.hartwig.hmftools.sage.read;

import static com.hartwig.hmftools.sage.read.ReadContextMatch.FULL;
import static com.hartwig.hmftools.sage.read.ReadContextMatch.NONE;
import static com.hartwig.hmftools.sage.read.ReadContextMatch.PARTIAL;

import java.util.Arrays;

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
    private final int repeatCount;
    private final String microhomology;

    public ReadContext(final String repeat, final int refPosition, final int readIndex, final int leftCentreIndex,
            final int rightCentreIndex, final int flankSize, final byte[] readBases) {
        assert (leftCentreIndex >= 0);
        assert (rightCentreIndex >= leftCentreIndex);

        this.position = refPosition;
        this.flankSize = flankSize;
        this.leftCentreIndex = leftCentreIndex;
        this.rightCentreIndex = rightCentreIndex;
        this.readBases = readBases;
        this.readIndex = readIndex;
        this.distance = 0;
        this.distanceCigar = Strings.EMPTY;
        this.repeat = repeat;
        this.microhomology = Strings.EMPTY;
        this.repeatCount = 0;
    }

    ReadContext(final String microhomology, int repeatCount, final String repeat, final int refPosition, final int readIndex,
            final int leftCentreIndex, final int rightCentreIndex, final int flankSize, byte[] refBases, @NotNull final SAMRecord record) {
        assert (leftCentreIndex >= 0);
        assert (rightCentreIndex >= leftCentreIndex);

        int recordLeftFlankStartIndex = Math.max(0, leftCentreIndex - flankSize);
        int recordLeftFlankLength = leftCentreIndex - recordLeftFlankStartIndex;

        int recordRightFlankEndIndex = Math.min(record.getReadBases().length - 1, rightCentreIndex + flankSize);
        int recordRightFlankLength = recordRightFlankEndIndex - rightCentreIndex;

        this.position = refPosition;
        this.flankSize = flankSize;
        this.leftCentreIndex = recordLeftFlankLength;
        this.rightCentreIndex = this.leftCentreIndex + rightCentreIndex - leftCentreIndex;
        this.readIndex = this.leftCentreIndex + readIndex - leftCentreIndex;
        this.repeat = repeat;
        this.repeatCount = repeatCount;
        this.microhomology = microhomology;
        this.readBases = Arrays.copyOfRange(record.getReadBases(), recordLeftFlankStartIndex, recordRightFlankEndIndex + 1);

        ReadContextDistance distance = new ReadContextDistance(recordLeftFlankStartIndex, recordRightFlankEndIndex, record, refBases);
        this.distance = distance.distance();
        this.distanceCigar = distance.cigar();
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
        return isComplete() && other.isComplete() && centreMatch(other.readIndex, other.readBases)
                && leftFlankMatchingBases(other.readIndex, other.readBases) == flankSize
                && rightFlankMatchingBases(other.readIndex, other.readBases) == flankSize;
    }

    int minCentreQuality(int readIndex, SAMRecord record) {
        int leftOffset = this.readIndex - leftCentreIndex;
        int rightOffset = rightCentreIndex - this.readIndex;

        int leftIndex = readIndex - leftOffset;
        int rightIndex = readIndex + rightOffset;

        int quality = Integer.MAX_VALUE;
        for (int i = leftIndex; i <= rightIndex; i++) {
            quality = Math.min(quality, record.getBaseQualities()[i]);
        }
        return quality;
    }

    public boolean phased(int offset, @NotNull final ReadContext other) {
        int otherReadIndex = other.readIndex + offset;

        boolean centreMatch = centreMatch(otherReadIndex, other.readBases);
        if (!centreMatch) {
            return false;
        }

        boolean otherCentreMatch = other.centreMatch(readIndex - offset, readBases);
        if (!otherCentreMatch) {
            return false;
        }

        int leftFlankingBases = leftFlankMatchingBases(otherReadIndex, other.readBases, flankSize);
        if (leftFlankingBases < 0) {
            return false;
        }

        int rightFlankingBases = rightFlankMatchingBases(otherReadIndex, other.readBases, flankSize);
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

        if (otherReadIndex < 0 || !isComplete()) {
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

    @VisibleForTesting
    boolean centreMatch(int otherRefIndex, byte[] otherBases) {

        int otherLeftCentreIndex = otherLeftCentreIndex(otherRefIndex);
        if (otherLeftCentreIndex < 0) {
            return false;
        }

        int otherRightCentreIndex = otherRightCentreIndex(otherRefIndex);
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

    private int otherLeftCentreIndex(int otherRefIndex) {
        return otherRefIndex + leftCentreIndex - readIndex;
    }

    private int otherRightCentreIndex(int otherRefIndex) {
        return otherRefIndex + rightCentreIndex - readIndex;
    }

    @VisibleForTesting
    int leftFlankMatchingBases(int otherRefIndex, byte[] otherBases) {
        return leftFlankMatchingBases(otherRefIndex, otherBases, flankSize);
    }

    private int leftFlankMatchingBases(int otherRefIndex, byte[] otherBases, int flankSize) {
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

    private int rightFlankMatchingBases(int otherRefIndex, byte[] otherBases, int flankSize) {
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

    public int leftFlankStartIndex() {
        return Math.max(0, leftCentreIndex - flankSize);
    }

    private int leftFlankLength() {
        return leftCentreIndex - leftFlankStartIndex();
    }

    public int rightFlankEndIndex() {
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

    public int repeatCount() {
        return repeatCount;
    }

    public byte[] readBases() {
        return readBases;
    }

    @NotNull
    public String mnvAdditionalAlt(int length) {
        return new String(readBases, readIndex - length + 1, length);
    }


}
