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
    private double jitterPenalty;
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
        return Math.max(0, quality - (int) jitterPenalty);
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
                MatchType match = readContext.matchAtPosition(readIndex, record.getReadBases());
                if (!match.equals(MatchType.NONE)) {
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
                    }
                } else {
                    final RealignedContext context = new Realigned().realigned(readContext, record.getReadBases());
                    switch (context.type()) {
                        case EXACT:
                            realigned++;
                            break;
                        case LENGTHENED:
                            jitterPenalty += jitterPenalty(context);
                            lengthened++;
                            break;
                        case SHORTENED:
                            jitterPenalty += jitterPenalty(context);
                            shortened++;
                            break;
                    }
                }
            }
        }
    }

    private void incrementQualityScores(int readBaseIndex, final SAMRecord record) {
        final int baseQuality = baseQuality(readBaseIndex, record);
        final int distanceFromReadEdge = readContext.distanceFromReadEdge(readBaseIndex, record);

        final int mapQuality = record.getMappingQuality();

        int modifiedMapQuality = mapQuality - 24 - 5 * (readContext.distance() - 1) - 15 * (record.getProperPairFlag() ? 0 : 1);
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
        return (int) jitterPenalty;
    }

    private double jitterPenalty(RealignedContext context) {
        return (0.25 * Math.max(0, context.repeatCount() - 3));
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
