package com.hartwig.hmftools.sage.read;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.sage.context.Realigned;
import com.hartwig.hmftools.sage.context.RealignedContext;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class ReadContextCounter implements GenomePosition, Consumer<SAMRecord> {
    private final VariantHotspot variant;
    private final ReadContext readContext;

    private int full;
    private int partial;
    private int realigned;
    private int lengthened;
    private int shortened;
    private int coverage;

    private int quality;
    private double jitterPenalty;
    private int baseQuality;
    private int mapQuality;

    private int excessiveInferredSize;
    private int inconsistentChromosome;
    private int improperPair;

    public ReadContextCounter(@NotNull final VariantHotspot variant, @NotNull final ReadContext readContext) {
        assert (readContext.isComplete());
        this.variant = variant;
        this.readContext = readContext;
    }

    @NotNull
    @Override
    public String chromosome() {
        return variant.chromosome();
    }

    @Override
    public long position() {
        return variant.position();
    }

    public int support() {
        return full + partial + realigned;
    }

    public int coverage() {
        return coverage;
    }

    public double vaf() {
        return coverage == 0 ? 0d : (double) support() / coverage;
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

        if (record.getAlignmentStart() <= variant.position() && record.getAlignmentEnd() >= variant.position()
                && readContext.isComplete()) {

            boolean covered = false;
            if (coverage < 1000) {
                int readIndex = record.getReadPositionAtReferencePosition(readContext.position()) - 1;
                if (readContext.isCentreCovered(readIndex, record.getReadBases())) {
                    covered = true;
                }

                ReadContextMatch match = readContext.matchAtPosition(readIndex, record.getReadBases());
                if (!match.equals(ReadContextMatch.NONE)) {
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
                            covered = true;
                            break;
                        case LENGTHENED:
                            jitterPenalty += jitterPenalty(context);
                            lengthened++;
                            covered = true;
                            break;
                        case SHORTENED:
                            jitterPenalty += jitterPenalty(context);
                            shortened++;
                            covered = true;
                            break;
                    }
                }
            }

            if (covered) {
                coverage++;
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
        return variant.ref().length() == variant.alt().length()
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
