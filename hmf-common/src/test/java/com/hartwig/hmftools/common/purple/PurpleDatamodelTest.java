package com.hartwig.hmftools.common.purple;

import static org.junit.Assert.assertNotNull;

import java.util.Optional;

import com.hartwig.hmftools.common.cobalt.ImmutableCobaltRatio;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.copynumber.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.sv.ImmutableStructuralVariantLegPloidy;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;
import com.hartwig.hmftools.common.purple.region.ImmutableEnrichedRegion;
import com.hartwig.hmftools.common.purple.region.ImmutableFittedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantImpl;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantLegImpl;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PurpleDatamodelTest {

    private static final String CHROMOSOME = "1";

    @Test
    public void testDefaultFittedRegion() {
        assertNotNull(createDefaultFittedRegion(CHROMOSOME, 1, 100).build());
    }

    @Test
    public void testDefaultCopyNumber() {
        assertNotNull(createCopyNumber(CHROMOSOME, 1, 100, 2).build());
    }

    @NotNull
    public static ImmutablePurpleCopyNumber.Builder createCopyNumber(@NotNull final String chromosome, final long start, final long end,
            final double copyNumber) {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome(chromosome)
                .start(start)
                .end(end)
                .averageTumorCopyNumber(copyNumber)
                .segmentStartSupport(SegmentSupport.NONE)
                .segmentEndSupport(SegmentSupport.NONE)
                .method(CopyNumberMethod.UNKNOWN)
                .bafCount(0)
                .averageObservedBAF(0.5)
                .averageActualBAF(0.5);
    }

    @NotNull
    public static ImmutableEnrichedRegion.Builder createObservedRegion(@NotNull final String chromosome, final long start, final long end) {
        return ImmutableEnrichedRegion.builder()
                .observedBAF(0.5)
                .bafCount(1)
                .chromosome(chromosome)
                .start(start)
                .end(end)
                .observedTumorRatio(1)
                .depthWindowCount(1)
                .observedNormalRatio(1)
                .ratioSupport(true)
                .svCluster(false)
                .status(GermlineStatus.DIPLOID)
                .gcContent(0.93)
                .support(SegmentSupport.NONE);
    }

    @NotNull
    public static ImmutableFittedRegion.Builder createDefaultFittedRegion(@NotNull final String chromosome, final long start,
            final long end) {
        final ObservedRegion observedRegion = createObservedRegion(chromosome, start, end).build();
        return ImmutableFittedRegion.builder()
                .from(observedRegion)
                .tumorCopyNumber(2)
                .tumorBAF(0.5)
                .fittedBAF(0)
                .fittedTumorCopyNumber(0)
                .deviation(0)
                .minorAllelePloidy(0)
                .minorAllelePloidyDeviation(0)
                .majorAllelePloidy(0)
                .majorAllelePloidyDeviation(0)
                .refNormalisedCopyNumber(2)
                .ratioSupport(true)
                .support(SegmentSupport.NONE)
                .ploidyPenalty(0);
    }

    @NotNull
    public static ImmutableStructuralVariantLegImpl.Builder createStartLeg(@NotNull final String startChromosome, final long startPosition,
            @NotNull final StructuralVariantType type) {

        final byte startOrientation;
        switch (type) {
            case DUP:
                startOrientation = -1;
                break;
            case BND:
            case INV:
                startOrientation = 1;
                break;
            default:
                startOrientation = 1;
                break;
        }

        return ImmutableStructuralVariantLegImpl.builder()
                .chromosome(startChromosome)
                .position(startPosition)
                .homology("")
                .orientation(startOrientation);
    }

    @NotNull
    public static ImmutableStructuralVariantLegImpl.Builder createEndLeg(@NotNull final String endChromosome, final long endPosition,
            @NotNull final StructuralVariantType type) {

        final byte endOrientation;
        switch (type) {
            case DUP:
                endOrientation = 1;
                break;
            case BND:
            case INV:
                endOrientation = 1;
                break;
            default:
                endOrientation = -1;
                break;
        }

        return ImmutableStructuralVariantLegImpl.builder()
                .chromosome(endChromosome)
                .position(endPosition)
                .orientation(endOrientation)
                .homology("");
    }

    @NotNull
    public static ImmutableStructuralVariantImpl.Builder createStructuralVariant(@NotNull final String startChromosome,
            final long startPosition, @NotNull final String endChromosome, final long endPosition,
            @NotNull final StructuralVariantType type, double startVaf, double endVaf) {

        return ImmutableStructuralVariantImpl.builder()
                .id("")
                .insertSequence("")
                .type(type)
                .start(createStartLeg(startChromosome, startPosition, type).alleleFrequency(startVaf).build())
                .end(createEndLeg(endChromosome, endPosition, type).alleleFrequency(endVaf).build())
                .imprecise(true)
                .somaticScore(0);
    }

    @NotNull
    public static ImmutableStructuralVariantImpl.Builder createStructuralVariant(@NotNull final String startChromosome,
            final long startPosition, @NotNull final String endChromosome, final long endPosition,
            @NotNull final StructuralVariantType type) {

        return ImmutableStructuralVariantImpl.builder()
                .id("")
                .insertSequence("")
                .type(type)
                .start(createStartLeg(startChromosome, startPosition, type).build())
                .end(createEndLeg(endChromosome, endPosition, type).build())
                .imprecise(true)
                .somaticScore(0);
    }
    @NotNull
    public static ImmutableStructuralVariantImpl.Builder createStructuralVariantSingleBreakend(@NotNull final String startChromosome,
           final long startPosition, double startVaf) {
        return ImmutableStructuralVariantImpl.builder()
                .id("")
                .insertSequence("")
                .type(StructuralVariantType.BND)
                .start(createStartLeg(startChromosome, startPosition, StructuralVariantType.BND).alleleFrequency(startVaf).build())
                .imprecise(false)
                .somaticScore(0);
    }

    @NotNull
    public static ImmutableCobaltRatio.Builder cobalt(@NotNull final String chromosome, long position, double ratio) {
        return ImmutableCobaltRatio.builder()
                .chromosome(chromosome)
                .position(position)
                .tumorReadCount(0)
                .referenceReadCount(0)
                .referenceGCRatio(1)
                .referenceGCDiploidRatio(1)
                .tumorGCRatio(ratio);
    }

    @NotNull
    public static ImmutableStructuralVariantLegPloidy.Builder svLegPloidy(int orientation, @NotNull final Optional<Double> leftCopyNumber,
            @NotNull final Optional<Double> rightCopyNumber, double ploidy) {
        return ImmutableStructuralVariantLegPloidy.builder()
                .chromosome(CHROMOSOME)
                .position(1)
                .orientation((byte) orientation)
                .vaf(0.5)
                .alleleFrequency(0.5)
                .homology("")
                .weight(1)
                .averageImpliedPloidy(ploidy)
                .unweightedImpliedPloidy(ploidy)
                .leftCopyNumber(leftCopyNumber)
                .rightCopyNumber(rightCopyNumber);
    }

}
