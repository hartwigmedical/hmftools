package com.hartwig.hmftools.common.purple;

import java.util.Collection;
import java.util.Random;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cobalt.ImmutableCobaltRatio;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.copynumber.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.purity.ImmutableFittedPurity;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;
import com.hartwig.hmftools.common.purple.region.ImmutableEnrichedRegion;
import com.hartwig.hmftools.common.purple.region.ImmutableFittedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.sv.ImmutableStructuralVariantImpl;
import com.hartwig.hmftools.common.sv.ImmutableStructuralVariantLegImpl;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class PurpleTestUtils
{
    private static final String CHROMOSOME = "1";

    public static double nextDouble(@NotNull final Random random) {
        return Math.round(random.nextDouble() * 10000D) / 10000D;
    }

    @NotNull
    public static ImmutableFittedPurity.Builder createRandomPurityBuilder(@NotNull Random random) {
        return ImmutableFittedPurity.builder()
                .purity(nextDouble(random))
                .normFactor(nextDouble(random))
                .score(nextDouble(random))
                .diploidProportion(nextDouble(random))
                .ploidy(nextDouble(random))
                .somaticPenalty(nextDouble(random));
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
                .depthWindowCount(1)
                .gcContent(0)
                .minStart(start)
                .maxStart(start)
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
                .unnormalisedObservedNormalRatio(1)
                .ratioSupport(true)
                .svCluster(false)
                .minStart(0)
                .maxStart(0)
                .status(GermlineStatus.DIPLOID)
                .gcContent(0.93)
                .support(SegmentSupport.NONE);
    }

    @NotNull
    public static ImmutableStructuralVariantImpl.Builder createStructuralVariant(@NotNull final String startChromosome,
            final long startPosition, @NotNull final String endChromosome, final long endPosition,
            @NotNull final StructuralVariantType type, double startVaf, double endVaf) {
        return ImmutableStructuralVariantImpl.builder()
                .id(Strings.EMPTY)
                .insertSequence(Strings.EMPTY)
                .insertSequenceAlignments(Strings.EMPTY)
                .type(type)
                .hotspot(false)
                .recovered(false)
                .qualityScore(0)
                .start(createStartLeg(startChromosome, startPosition, type).alleleFrequency(startVaf).build())
                .end(createEndLeg(endChromosome, endPosition, type).alleleFrequency(endVaf).build())
                .startContext(dummyContext())
                .imprecise(true);
    }

    @NotNull
    public static ImmutableStructuralVariantImpl.Builder createStructuralVariant(@NotNull final String startChromosome,
            final long startPosition, @NotNull final String endChromosome, final long endPosition,
            @NotNull final StructuralVariantType type) {
        return ImmutableStructuralVariantImpl.builder()
                .id(Strings.EMPTY)
                .insertSequence(Strings.EMPTY)
                .insertSequenceAlignments(Strings.EMPTY)
                .type(type)
                .qualityScore(0)
                .hotspot(false)
                .recovered(false)
                .start(createStartLeg(startChromosome, startPosition, type).build())
                .end(createEndLeg(endChromosome, endPosition, type).build())
                .startContext(dummyContext() )
                .imprecise(true);
    }

    @NotNull
    public static ImmutableStructuralVariantImpl.Builder createStructuralVariantSingleBreakend(@NotNull final String startChromosome,
            final long startPosition, double startVaf) {
        return ImmutableStructuralVariantImpl.builder()
                .id(Strings.EMPTY)
                .insertSequence(Strings.EMPTY)
                .insertSequenceAlignments(Strings.EMPTY)
                .qualityScore(0)
                .hotspot(false)
                .recovered(false)
                .type(StructuralVariantType.BND)
                .start(createStartLeg(startChromosome, startPosition, StructuralVariantType.BND).alleleFrequency(startVaf).build())
                .startContext(dummyContext())
                .imprecise(false);
    }

    @NotNull
    public static ImmutableFittedRegion.Builder createDefaultFittedRegion(@NotNull final String chromosome, final long start,
            final long end) {
        final ObservedRegion observedRegion = PurpleTestUtils.createObservedRegion(chromosome, start, end).build();
        return ImmutableFittedRegion.builder()
                .from(observedRegion)
                .tumorCopyNumber(2)
                .tumorBAF(0.5)
                .fittedBAF(0)
                .fittedTumorCopyNumber(0)
                .deviationPenalty(0)
                .minorAlleleCopyNumberDeviation(0)
                .majorAlleleCopyNumberDeviation(0)
                .refNormalisedCopyNumber(2)
                .ratioSupport(true)
                .support(SegmentSupport.NONE)
                .eventPenalty(0);
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
                .homology(Strings.EMPTY)
                .orientation(startOrientation)
                .anchoringSupportDistance(0);
    }

    @NotNull
    private static ImmutableStructuralVariantLegImpl.Builder createEndLeg(@NotNull final String endChromosome, final long endPosition,
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
                .homology("")
                .anchoringSupportDistance(0);
    }

    @NotNull
    private static VariantContext dummyContext() {
        final Collection<Allele> alleles = Lists.newArrayList(Allele.create("N", true));

        return new VariantContextBuilder("purple", "2", 1, 1, alleles)
                .noGenotypes()
                .make();
    }

}
