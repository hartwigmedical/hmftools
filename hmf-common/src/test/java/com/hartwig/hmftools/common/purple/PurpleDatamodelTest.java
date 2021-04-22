package com.hartwig.hmftools.common.purple;

import static org.junit.Assert.assertNotNull;

import java.util.Collection;
import java.util.Optional;

import com.google.common.collect.Lists;
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

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class PurpleDatamodelTest {

    private static final String CHROMOSOME = "1";

    @Test
    public void testDefaultFittedRegion() {
        assertNotNull(createDefaultFittedRegion(CHROMOSOME, 1, 100).build());
    }

    @Test
    public void testDefaultCopyNumber() {
        assertNotNull(PurpleTestUtils.createCopyNumber(CHROMOSOME, 1, 100, 2).build());
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
    public static ImmutableStructuralVariantLegPloidy.Builder svLegPloidy(int orientation, @NotNull final Optional<Double> leftCopyNumber,
            @NotNull final Optional<Double> rightCopyNumber, double ploidy) {
        return ImmutableStructuralVariantLegPloidy.builder()
                .chromosome(CHROMOSOME)
                .position(1)
                .orientation((byte) orientation)
                .observedVaf(0.5)
                .adjustedVaf(0.5)
                .alleleFrequency(0.5)
                .homology("")
                .anchoringSupportDistance(0)
                .weight(1)
                .averageImpliedPloidy(ploidy)
                .unweightedImpliedPloidy(ploidy)
                .leftCopyNumber(leftCopyNumber)
                .rightCopyNumber(rightCopyNumber);
    }

    @NotNull
    private static VariantContext dummyContext() {
        final Collection<Allele> alleles = Lists.newArrayList(Allele.create("N", true));

        return new VariantContextBuilder("purple", "2", 1, 1, alleles)
                .noGenotypes()
                .make();
    }
}
