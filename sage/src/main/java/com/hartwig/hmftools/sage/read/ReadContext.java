package com.hartwig.hmftools.sage.read;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class ReadContext {

    private final int position;
    private final int distance;
    private final String distanceCigar;
    private final String repeat;
    private final int repeatCount;
    private final String microhomology;
    private final IndexedBases readBases;
    private final IndexedBases refBases;

    @VisibleForTesting
    ReadContext(final String repeat, final int refPosition, final int readIndex, final int leftCentreIndex, final int rightCentreIndex,
            final int flankSize, final byte[] readBases, final String microhomology) {
        assert (leftCentreIndex >= 0);
        assert (rightCentreIndex >= leftCentreIndex);

        this.position = refPosition;
        this.distance = 0;
        this.distanceCigar = Strings.EMPTY;
        this.repeat = repeat;
        this.microhomology = microhomology;
        this.repeatCount = 0;
        this.readBases = new IndexedBases(refPosition, readIndex, leftCentreIndex, rightCentreIndex, flankSize, readBases);
        this.refBases = new IndexedBases(refPosition, readIndex, leftCentreIndex, rightCentreIndex, flankSize, readBases);
    }

    ReadContext(final String microhomology, int repeatCount, final String repeat, final int refPosition, final int readIndex,
            final int leftCentreIndex, final int rightCentreIndex, final int flankSize, @NotNull final IndexedBases refSequence,
            @NotNull final SAMRecord record) {
        assert (leftCentreIndex >= 0);
        assert (rightCentreIndex >= leftCentreIndex);

        this.position = refPosition;
        this.repeat = repeat;
        this.repeatCount = repeatCount;
        this.microhomology = microhomology;

        int recordLeftFlankStartIndex = Math.max(0, leftCentreIndex - flankSize);
        int recordRightFlankEndIndex = Math.min(record.getReadBases().length - 1, rightCentreIndex + flankSize);

        ReadContextDistance distance = new ReadContextDistance(recordLeftFlankStartIndex, recordRightFlankEndIndex, record, refSequence);
        this.distance = distance.distance();
        this.distanceCigar = distance.cigar();

        this.readBases = new IndexedBases(refPosition, readIndex, leftCentreIndex, rightCentreIndex, flankSize, record.getReadBases());

        int refIndex = refSequence.index(position);
        this.refBases = IndexedBases.resize(position,
                refIndex,
                refIndex + leftCentreIndex - readIndex,
                refIndex + rightCentreIndex - readIndex,
                0,
                refSequence.bases());
    }

    private ReadContext(@NotNull final ReadContext clone) {
        this.position = clone.position;
        this.repeat = clone.repeat;
        this.repeatCount = clone.repeatCount;
        this.microhomology = clone.microhomology;
        this.distance = clone.distance;
        this.distanceCigar = clone.distanceCigar;
        this.refBases = clone.refBases;

        this.readBases = IndexedBases.resize(position,
                clone.readBases.index(),
                clone.readBases.leftCentreIndex(),
                clone.readBases.rightCentreIndex(),
                clone.readBases.flankSize(),
                clone.readBases.bases());

    }

    @NotNull
    public ReadContext minimiseFootprint() {
        return new ReadContext(this);
    }

    public int position() {
        return position;
    }

    public boolean isComplete() {
        return readBases.flanksComplete();
    }

    public boolean isFullMatch(@NotNull final ReadContext other) {
        return isComplete() && other.isComplete() && readBases.coreMatch(other.readIndex(), other.readBases())
                && readBases.leftFlankMatchingBases(other.readIndex(), other.readBases()) == flankSize()
                && readBases.rightFlankMatchingBases(other.readIndex(), other.readBases()) == flankSize();
    }

    int minCentreQuality(int readIndex, SAMRecord record) {
        int leftOffset = this.readIndex() - readBases.leftCentreIndex();
        int rightOffset = readBases.rightCentreIndex() - this.readIndex();

        int leftIndex = readIndex - leftOffset;
        int rightIndex = readIndex + rightOffset;

        int quality = Integer.MAX_VALUE;
        for (int i = leftIndex; i <= rightIndex; i++) {
            quality = Math.min(quality, record.getBaseQualities()[i]);
        }
        return quality;
    }

    int avgCentreQuality(int readIndex, @NotNull final SAMRecord record) {
        int leftOffset = this.readIndex() - readBases.leftCentreIndex();
        int rightOffset = readBases.rightCentreIndex() - this.readIndex();

        int leftIndex = readIndex - leftOffset;
        int rightIndex = readIndex + rightOffset;

        float quality = 0;
        for (int i = leftIndex; i <= rightIndex; i++) {
            quality += record.getBaseQualities()[i];
        }
        return Math.round(quality / (rightIndex - leftIndex + 1));
    }

    public boolean phased(int offset, @NotNull final ReadContext other) {
        return readBases.phased(offset, other.readBases);
    }

    public boolean isCentreCovered(int otherReadIndex, byte[] otherBases) {
        return readBases.isCentreCovered(otherReadIndex, otherBases);
    }

    @NotNull
    public ReadContextMatch matchAtPosition(int otherReadIndex, byte[] otherBases) {
        return readBases.matchAtPosition(otherReadIndex, otherBases);
    }

    public int readBasesPositionIndex() {
        return readBases.index();
    }

    public int readBasesLeftFlankIndex() {
        return readBases.leftFlankIndex();
    }

    public int readBasesRightFlankIndex() {
        return readBases.rightFlankIndex();
    }

    public int readBasesLeftCentreIndex() {
        return readBases.leftCentreIndex();
    }

    public int readBasesRightCentreIndex() {
        return readBases.rightCentreIndex();
    }

    @Override
    public String toString() {
        return readBases.centerString();
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
        return readBases.centerString();
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
        return readBases.bases();
    }

    private int readIndex() {
        return readBases.index();
    }

    private int flankSize() {
        return readBases.flankSize();
    }

    public byte[] refTrinucleotideContext(int position) {
        return refBases.trinucleotideContext(position);
    }

}
