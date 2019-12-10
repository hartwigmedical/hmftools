package com.hartwig.hmftools.sage.read;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.QualityConfig;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.context.Realigned;
import com.hartwig.hmftools.sage.context.RealignedContext;
import com.hartwig.hmftools.sage.context.RealignedType;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class ReadContextCounter implements GenomePosition {
    private final VariantHotspot variant;
    private final ReadContext readContext;

    private int full;
    private int partial;
    private int core;
    private int reference;
    private int realigned;
    private int lengthened;
    private int shortened;
    private int coverage;

    private int quality;

    private double jitterPenalty;

    private int improperPair;

    public ReadContextCounter(@NotNull final VariantHotspot variant, @NotNull final ReadContext readContext) {
        assert (readContext.isComplete());
        this.variant = variant;
        this.readContext = readContext;
    }

    public ReadContextCounter(@NotNull final VariantHotspot variant, @NotNull final ReadContextCounter left,
            @NotNull final ReadContextCounter right) {
        this.variant = variant;
        this.readContext = left.readContext;
        this.full = Math.min(left.full, right.full);
        this.partial = Math.min(left.partial, right.partial);
        this.core = Math.min(left.core, right.core);
        this.reference = Math.min(left.reference, right.reference);
        this.realigned = Math.min(left.realigned, right.realigned);
        this.lengthened = Math.min(left.lengthened, right.lengthened);
        this.shortened = Math.min(left.shortened, right.shortened);
        this.coverage = Math.min(left.coverage, right.coverage);
        this.quality = Math.min(left.quality, right.quality);
        this.jitterPenalty = Math.min(left.jitterPenalty, right.jitterPenalty);
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

    public int altSupport() {
        return full + partial + core + realigned;
    }

    public int coverage() {
        return coverage;
    }

    public int refSupport() {
        return reference;
    }

    public int depth() {
        return coverage;
    }

    public double af(double support) {
        return coverage == 0 ? 0d : support / depth();
    }

    public double vaf() {
        return af(altSupport());
    }

    public int quality() {
        return Math.max(0, quality - (int) jitterPenalty);
    }

    public int[] rcc() {
        return new int[] { full, partial, core, realigned, shortened, lengthened, reference, coverage };
    }

    public int[] qual() {
        return new int[] { quality(), 0, 0, qualityJitterPenalty() };
    }

    public int improperPair() {
        return improperPair;
    }

    @NotNull
    public ReadContext readContext() {
        return readContext;
    }

    @Override
    public String toString() {
        return readContext.toString();
    }

    public void accept(final boolean realign, final SAMRecord record, final SageConfig sageConfig, IndexedBases refSequence) {
        final QualityConfig qualityConfig = sageConfig.qualityConfig();

        try {

            if (record.getAlignmentStart() <= variant.position() && record.getAlignmentEnd() >= variant.position()
                    && readContext.isComplete()) {

                if (coverage >= sageConfig.maxReadDepth()) {
                    return;
                }

                int readIndex = record.getReadPositionAtReferencePosition(readContext.position()) - 1;

                // Check if FULL, PARTIAL, OR CORE
                final ReadContextMatch match = readContext.matchAtPosition(readIndex, record.getReadBases());
                if (!match.equals(ReadContextMatch.NONE)) {
                    switch (match) {
                        case FULL:
                            full++;
                            incrementQualityFlags(record);
                            incrementQualityScores(readIndex, record, qualityConfig);
                            break;
                        case PARTIAL:
                            partial++;
                            incrementQualityFlags(record);
                            incrementQualityScores(readIndex, record, qualityConfig);
                            break;
                        case CORE:
                            core++;
                            break;
                    }

                    coverage++;
                    return;
                }

                // Check if realigned
                final RealignedContext realignment = realignmentContext(realign, readIndex, record);
                if (realignment.type().equals(RealignedType.EXACT)) {
                    realigned++;
                    coverage++;
                    return;
                }

                // Check if lengthened, shortened AND/OR reference!
                boolean covered = readContext.isCentreCovered(readIndex, record.getReadBases());

                switch (realignment.type()) {
                    case LENGTHENED:
                        jitterPenalty += qualityConfig.jitterPenalty(realignment.repeatCount());
                        lengthened++;
                        covered = true;
                        break;
                    case SHORTENED:
                        jitterPenalty += qualityConfig.jitterPenalty(realignment.repeatCount());
                        shortened++;
                        covered = true;
                        break;
                }

                if (covered && readContext.matchesRef(readIndex, record.getReadBases())) {
                    reference++;
                }

                if (covered) {
                    coverage++;
                }
            }
        } catch (Exception e) {
            System.out.println("Error at chromosome: " + chromosome() + ", position: " + position());
            throw e;
        }
    }

    @NotNull
    private RealignedContext realignmentContext(boolean realign, int readIndex, SAMRecord record) {
        if (!realign) {
            return new RealignedContext(RealignedType.NONE, 0);
        }

        if (variant.isSNV() && record.getCigar().getCigarElements().size() == 1) {
            return new Realigned().realignedAroundIndex(readContext, readIndex, record.getReadBases());
        }

        return new Realigned().realignedInEntireRecord(readContext, record.getReadBases());
    }

    private void incrementQualityScores(int readBaseIndex, final SAMRecord record, final QualityConfig qualityConfig) {
        final int baseQuality = baseQuality(readBaseIndex, record);
        final int distanceFromReadEdge = readContext.distanceFromReadEdge(readBaseIndex, record);

        final int mapQuality = record.getMappingQuality();

        int modifiedMapQuality = qualityConfig.modifiedMapQuality(mapQuality, readContext.distance(), record.getProperPairFlag());
        int modifiedBaseQuality = qualityConfig.modifiedBaseQuality(baseQuality, distanceFromReadEdge);

        this.quality += Math.max(0, Math.min(modifiedMapQuality, modifiedBaseQuality));
    }

    private int baseQuality(int readBaseIndex, SAMRecord record) {
        return variant.ref().length() == variant.alt().length()
                ? record.getBaseQualities()[readBaseIndex]
                : readContext.minCentreQuality(readBaseIndex, record);
    }

    private int qualityJitterPenalty() {
        return (int) jitterPenalty;
    }

    public boolean incrementCounters(@NotNull final ReadContext other) {
        if (readContext.isFullMatch(other)) {
            full++;
            return true;
        }

        return false;
    }

    private void incrementQualityFlags(@NotNull final SAMRecord record) {
        if (!record.getProperPairFlag()) {
            improperPair++;
        }
    }

    private int indelLength(@NotNull final SAMRecord record) {
        int result = 0;
        for (CigarElement cigarElement : record.getCigar()) {
            switch (cigarElement.getOperator()) {
                case I:
                case D:
                    result += cigarElement.getLength();
            }

        }

        return result;
    }

}
