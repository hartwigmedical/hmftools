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

    private final IndexedBases readIndexedBases;

//    private final int refIndex;
//    private final byte refBases[];

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
        this.readIndexedBases = new IndexedBases(readIndex, leftCentreIndex, rightCentreIndex, flankSize, readBases);
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
        this.readIndexedBases = new IndexedBases(this.readIndex, this.leftCentreIndex, this.rightCentreIndex, flankSize, readBases);
    }

    public int position() {
        return position;
    }

    int distanceFromReadEdge(int readIndex, SAMRecord record) {
        int leftOffset = this.readIndex() - leftCentreIndex;
        int rightOffset = rightCentreIndex - this.readIndex();

        int leftIndex = readIndex - leftOffset;
        int rightIndex = readIndex + rightOffset;

        return Math.min(leftIndex, record.getReadBases().length - rightIndex - 1);
    }

    public boolean isComplete() {
        return readIndexedBases.flanksComplete();
    }

    public boolean isFullMatch(@NotNull final ReadContext other) {
        return isComplete() && other.isComplete() && centreMatch(other.readIndex, other.readBases())
                && leftFlankMatchingBases(other.readIndex, other.readBases()) == flankSize
                && rightFlankMatchingBases(other.readIndex, other.readBases()) == flankSize;
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
        return readIndexedBases.phased(offset, other.readIndexedBases);
    }

    public boolean isCentreCovered(int otherReadIndex, byte[] otherBases) {
        return readIndexedBases.isCentreCovered(otherReadIndex, otherBases);
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
        return readIndexedBases.centreMatch(otherRefIndex, otherBases);
    }

    private int otherLeftCentreIndex(int otherRefIndex) {
        return otherRefIndex + leftCentreIndex - readIndex;
    }

    private int otherRightCentreIndex(int otherRefIndex) {
        return otherRefIndex + rightCentreIndex - readIndex;
    }

    @VisibleForTesting
    int leftFlankMatchingBases(int otherRefIndex, byte[] otherBases) {
        return readIndexedBases.leftFlankMatchingBases(otherRefIndex, otherBases, flankSize);
    }

    @VisibleForTesting
    int rightFlankMatchingBases(int otherRefIndex, byte[] otherBases) {
        return readIndexedBases.rightFlankMatchingBases(otherRefIndex, otherBases, flankSize);
    }


    public int leftFlankStartIndex() {
        return Math.max(0, leftCentreIndex - flankSize);
    }


    public int rightFlankEndIndex() {
        return Math.min(readBases.length - 1, rightCentreIndex + flankSize);
    }


    @Override
    public String toString() {
        return readIndexedBases.toString();
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
        return readIndexedBases.centerString();
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
        return readIndexedBases.bases();
    }

    private int readIndex() {
        return readIndexedBases.index();
    }

    @NotNull
    public String mnvAdditionalAlt(int length) {
        return new String(readBases, readIndex - length + 1, length);
    }


}
