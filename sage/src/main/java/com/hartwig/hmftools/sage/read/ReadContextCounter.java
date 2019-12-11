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

    private int fullQuality;
    private int partialQuality;
    private int coreQuality;
    private int realignedQuality;
    private int referenceQuality;
    private int totalQuality;

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
        this.fullQuality = Math.min(left.fullQuality, right.fullQuality);
        this.partialQuality = Math.min(left.partialQuality, right.partialQuality);
        this.coreQuality = Math.min(left.coreQuality, right.coreQuality);
        this.realignedQuality = Math.min(left.realignedQuality, right.realignedQuality);
        this.referenceQuality = Math.min(left.referenceQuality, right.referenceQuality);
        this.totalQuality = Math.min(left.totalQuality, right.totalQuality);
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

    public int refSupport() {
        return reference;
    }

    public int coverage() {
        return coverage;
    }

    public int depth() {
        return coverage;
    }

    public double vaf() {
        return af(tumorQuality());
    }

    public double refAllelicFrequency() {
        return af(refQuality());
    }

    public int refQuality() {
        return referenceQuality;
    }

    private double af(double support) {
        return coverage == 0 ? 0d : support / totalQuality;
    }

    public int tumorQuality() {
        int tumorQuality = fullQuality + partialQuality + coreQuality + realignedQuality;
        return Math.max(0, tumorQuality - (int) jitterPenalty);
    }

    public int[] counts() {
        return new int[] { full, partial, core, realigned, reference, coverage };
    }

    public int[] jitter() {
        return new int[] { shortened, lengthened, qualityJitterPenalty() };
    }

    public int[] quality() {
        return new int[] { fullQuality, partialQuality, coreQuality, realignedQuality, referenceQuality, totalQuality };
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
                boolean covered = readContext.isCentreCovered(readIndex, record.getReadBases());
                if (!covered) {
                    return;
                }

                // TODO: Check if this is okay? Should we check for jitter if quality 0
                double quality = calculateQualityScore(readIndex, record, qualityConfig, refSequence);
                if (quality <= 0) {
                    return;
                }

                coverage++;
                totalQuality += quality;

                // Check if FULL, PARTIAL, OR CORE
                final ReadContextMatch match = readContext.matchAtPosition(readIndex, record.getReadBases());
                if (!match.equals(ReadContextMatch.NONE)) {
                    switch (match) {
                        case FULL:
                            incrementQualityFlags(record);
                            full++;
                            fullQuality += quality;
                            break;
                        case PARTIAL:
                            incrementQualityFlags(record);
                            partial++;
                            partialQuality += quality;
                            break;
                        case CORE:
                            incrementQualityFlags(record);
                            core++;
                            coreQuality += quality;
                            break;
                    }

                    return;
                }

                // Check if realigned
                final RealignedContext realignment = realignmentContext(realign, readIndex, record);
                if (realignment.type().equals(RealignedType.EXACT)) {
                    realigned++;
                    realignedQuality += quality;
                    return;
                }

                // Check if lengthened, shortened AND/OR reference!

                switch (realignment.type()) {
                    case LENGTHENED:
                        jitterPenalty += qualityConfig.jitterPenalty(realignment.repeatCount());
                        lengthened++;
                        break;
                    case SHORTENED:
                        jitterPenalty += qualityConfig.jitterPenalty(realignment.repeatCount());
                        shortened++;
                        break;
                }


                if (readContext.matchesRef(readIndex, record.getReadBases())) {
                    reference++;
                    referenceQuality += quality;
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

        int indelLength = indelLength(record);
        return new Realigned().realignedAroundIndex(readContext,
                readIndex,
                record.getReadBases(),
                Math.max(indelLength, Realigned.MAX_REPEAT_SIZE));
    }

    private double calculateQualityScore(int readBaseIndex, final SAMRecord record, final QualityConfig qualityConfig,
            final IndexedBases refSequence) {
        final int distanceFromRef = readDistanceFromRef(readBaseIndex, record, refSequence);

        final int baseQuality = baseQuality(readBaseIndex, record);
        final int distanceFromReadEdge = readDistanceFromEdge(readBaseIndex, record);

        final int mapQuality = record.getMappingQuality();

        int modifiedMapQuality = qualityConfig.modifiedMapQuality(mapQuality, distanceFromRef, record.getProperPairFlag());
        int modifiedBaseQuality = qualityConfig.modifiedBaseQuality(baseQuality, distanceFromReadEdge);

        return Math.max(0, Math.min(modifiedMapQuality, modifiedBaseQuality));
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

    private int readDistanceFromEdge(int readIndex, @NotNull final SAMRecord record) {
        int index = readContext.readBasesPositionIndex();
        int leftIndex = readContext.readBasesLeftCentreIndex();
        int rightIndex = readContext.readBasesRightCentreIndex();

        int leftOffset = index - leftIndex;
        int rightOffset = rightIndex - index;

        int adjustedLeftIndex = readIndex - leftOffset;
        int adjustedRightIndex = readIndex + rightOffset;

        return Math.max(0, Math.min(adjustedLeftIndex, record.getReadBases().length - 1 - adjustedRightIndex));
    }

    private int readDistanceFromRef(int readIndex, final SAMRecord record, final IndexedBases refSequence) {
        int index = readContext.readBasesPositionIndex();
        int leftIndex = readContext.readBasesLeftFlankIndex();
        int rightIndex = readContext.readBasesRightFlankIndex();

        int leftOffset = index - leftIndex;
        int rightOffset = rightIndex - index;

        return new ReadContextDistance(Math.max(0, readIndex - leftOffset),
                Math.min(record.getReadBases().length - 1, readIndex + rightOffset),
                record,
                refSequence).distance();
    }

}
