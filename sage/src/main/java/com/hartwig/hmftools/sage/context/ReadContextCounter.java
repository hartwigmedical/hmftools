package com.hartwig.hmftools.sage.context;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.position.GenomePosition;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class ReadContextCounter implements GenomePosition, Consumer<SAMRecord> {
    private final VariantHotspot hotspot;
    private final ReadContextImproved readContext;

    private int full;
    private int partial;
    private int realigned;

    private int quality;
    private int baseQuality;
    private int mapQuality;

    private int excessiveInferredSize;
    private int inconsistentChromosome;
    private int improperPair;

    private int coverage;

    public ReadContextCounter(@NotNull final VariantHotspot hotspot, @NotNull final ReadContextImproved readContext) {
        assert (readContext.isComplete());
        this.hotspot = hotspot;
        this.readContext = readContext;
    }

    @NotNull
    @Override
    public String chromosome() {
        return hotspot.chromosome();
    }

    @Override
    public long position() {
        return hotspot.position();
    }

    public int full() {
        return full;
    }

    public int partial() {
        return partial;
    }

    public int realigned() {
        return realigned;
    }

    public int quality() {
        return quality;
    }

    public int baseQuality() {
        return baseQuality;
    }

    public int mapQuality() {
        return mapQuality;
    }

    public int[] rcq() {
        return new int[] { improperPair, inconsistentChromosome, excessiveInferredSize };
    }

    public ReadContextImproved readContext() {
        return readContext;
    }

    @Override
    public String toString() {
        return readContext.toString();
    }

    @Override
    public void accept(final SAMRecord record) {

        if (record.getAlignmentStart() <= hotspot.position() && record.getAlignmentEnd() >= hotspot.position()
                && readContext.isComplete()) {
            coverage++;
            if (coverage < 1000) {
                int readIndex = record.getReadPositionAtReferencePosition(readContext.position()) - 1;
                if (readIndex >= 0 && incrementCounters(readIndex, record.getReadBases())) {
                    incrementQualityScores(readIndex, record);

                    if (Math.abs(record.getInferredInsertSize()) >= 1000) {
                        excessiveInferredSize++;
                    }

                    if (!record.getReferenceName().equals(record.getMateReferenceName())) {
                        inconsistentChromosome++;
                    }

                    if (!record.getProperPairFlag()) {
                        improperPair++;
                    }
                } else if (readContext.isWithin(record.getReadBases())) {
                    realigned++;
                }
            }
        }
    }

    private void incrementQualityScores(int readBasePosition, final SAMRecord record) {
        final int distanceFromReadEdge = Math.min(readBasePosition, record.getReadBases().length - readBasePosition - 1);
        final int mapQuality = record.getMappingQuality();
        final int baseQuality = record.getBaseQualities()[readBasePosition];
        this.mapQuality += mapQuality;
        this.baseQuality += baseQuality;
        this.quality += quality(mapQuality, baseQuality, distanceFromReadEdge);
    }

    private double quality(int mapQuality, int baseQuality, int distanceFromEdge) {
        final int quality = Math.min(Math.min(Math.max(0, mapQuality - 12), baseQuality), distanceFromEdge);
        return Math.max(0, quality - 12);
    }

    private boolean incrementCounters(int otherRefIndex, byte[] otherBases) {
        final ReadContextMatch match = readContext.matchAtPosition(otherRefIndex, otherBases);
        if (match == ReadContextMatch.FULL) {
            full++;
            return true;
        }

        if (match == ReadContextMatch.PARTIAL) {
            partial++;
            return true;
        }

        return false;
    }

    public boolean incrementCounters(@NotNull final ReadContextImproved other) {
        if (readContext.isFullMatch(other)) {
            full++;
            return true;
        }

        return false;
    }

}
