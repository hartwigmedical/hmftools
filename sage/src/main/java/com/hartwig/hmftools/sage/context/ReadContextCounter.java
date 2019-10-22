package com.hartwig.hmftools.sage.context;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.position.GenomePosition;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class ReadContextCounter implements GenomePosition, Consumer<SAMRecord> {
    private final VariantHotspot hotspot;
    private final ReadContext readContext;

    private int full;
    private int partial;
    private int realigned;
    private int lengthened;
    private int shortened;

    private int quality;
    private int baseQuality;
    private int mapQuality;

    private int excessiveInferredSize;
    private int inconsistentChromosome;
    private int improperPair;

    private int coverage;

    public ReadContextCounter(@NotNull final VariantHotspot hotspot, @NotNull final ReadContext readContext) {
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
        return quality - qualityJitterPenalty();
    }

    public int baseQuality() {
        return baseQuality;
    }

    public int mapQuality() {
        return mapQuality;
    }

    public int[] rcc() {
        return new int[] { full, partial, realigned, shortened, lengthened };
    }

    public int[] rcq() {
        return new int[] { improperPair, inconsistentChromosome, excessiveInferredSize };
    }

    @NotNull
    public ReadContext readContext() {
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
                ReadContextMatch match = readContext.matchAtPosition(readIndex, record.getReadBases());
                switch (match) {
                    case FULL:
                        full++;
                        incrementQualityFlags(record);
                        incrementQualityScores(readIndex, record);
                        break;
                    case PARTIAL:
                        partial++;
                        incrementQualityFlags(record);
                        incrementQualityScores(readIndex, record);
                        break;
                    case JITTER_REMOVED:
                        shortened++;
                        break;
                    case JITTER_ADDED:
                        lengthened++;
                        break;
                    default:
                        if (readContext.isWithin(record.getReadBases())) {
                            realigned++;
                        }
                }
            }
        }
    }

    private void incrementQualityScores(int readBaseIndex, final SAMRecord record) {
        final int baseQuality = baseQuality(readBaseIndex, record);
        final int distanceFromReadEdge = readContext.distanceFromReadEdge(readBaseIndex, record);

        final int mapQuality = record.getMappingQuality();

        int modifiedMapQuality = mapQuality - 24 - 5 * (readContext.distance() - 1) - 15 * (record.getProperPairFlag() ? 1 : 0);
        int modifiedBaseQuality = Math.min(baseQuality, distanceFromReadEdge) - 12;

        this.mapQuality += mapQuality;
        this.baseQuality += baseQuality;
        this.quality += Math.max(0, Math.min(modifiedMapQuality, modifiedBaseQuality));
    }

    private int baseQuality(int readBaseIndex, SAMRecord record) {
        return hotspot.ref().length() == hotspot.alt().length()
                ? record.getBaseQualities()[readBaseIndex]
                : readContext.minCentreQuality(readBaseIndex, record);
    }

    public int qualityJitterPenalty() {
        return (int) (0.25 * Math.max(0, readContext.repeatCount() - 3) * (lengthened + shortened));
    }

    public boolean incrementCounters(@NotNull final ReadContext other) {
        if (readContext.isFullMatch(other)) {
            full++;
            return true;
        }

        return false;
    }

    private void incrementQualityFlags(@NotNull final SAMRecord record) {
        if (Math.abs(record.getInferredInsertSize()) >= 1000) {
            excessiveInferredSize++;
        }

        if (!record.getReferenceName().equals(record.getMateReferenceName())) {
            inconsistentChromosome++;
        }

        if (!record.getProperPairFlag()) {
            improperPair++;
        }
    }

}
