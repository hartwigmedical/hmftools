package com.hartwig.hmftools.common.purple;

import java.util.Collection;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
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

    public static CobaltRatio cobalt(@NotNull final String chromosome, int position, double ratio)
    {
        return new CobaltRatio(
                chromosome,
                position,
                0, // referenceReadDepth
                1, // referenceGCRatio
                0.5, // referenceGcContent
                1, // referenceGCDiploidRatio
                0, // tumorReadDepth
                ratio, // tumorGCRatio
                0.5); // tumorGcContent
    }

    public static ImmutablePurpleCopyNumber.Builder createCopyNumber(
            final String chromosome, final int start, final int end, final double copyNumber,
            final SegmentSupport startSupport, final SegmentSupport endSupport)
    {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome(chromosome)
                .start(start)
                .end(end)
                .averageTumorCopyNumber(copyNumber)
                .segmentStartSupport(startSupport)
                .segmentEndSupport(endSupport)
                .method(CopyNumberMethod.UNKNOWN)
                .bafCount(0)
                .depthWindowCount(1)
                .gcContent(0)
                .minStart(start)
                .maxStart(start)
                .averageObservedBAF(0.5)
                .averageActualBAF(0.5);
    }

    public static ImmutablePurpleCopyNumber.Builder createCopyNumber(
            final String chromosome, final int start, final int end, final double copyNumber)
    {
        return createCopyNumber(chromosome, start, end, copyNumber, SegmentSupport.NONE, SegmentSupport.NONE);
    }

    public static ImmutableStructuralVariantImpl.Builder createStructuralVariant(
            final String startChromosome, final int startPosition, final String endChromosome, final int endPosition,
            final StructuralVariantType type, double startVaf, double endVaf)
    {
        return ImmutableStructuralVariantImpl.builder()
                .id(Strings.EMPTY)
                .insertSequence(Strings.EMPTY)
                .insertSequenceAlignments(Strings.EMPTY)
                .type(type)
                .hotspot(false)
                .qualityScore(0)
                .start(createStartLeg(startChromosome, startPosition, type).alleleFrequency(startVaf).build())
                .end(createEndLeg(endChromosome, endPosition, type).alleleFrequency(endVaf).build())
                .startContext(dummyContext());
    }

    public static ImmutableStructuralVariantImpl.Builder createStructuralVariant(
            final String startChromosome, final int startPosition, final String endChromosome, final int endPosition,
            final StructuralVariantType type)
    {
        return ImmutableStructuralVariantImpl.builder()
                .id(Strings.EMPTY)
                .insertSequence(Strings.EMPTY)
                .insertSequenceAlignments(Strings.EMPTY)
                .type(type)
                .qualityScore(0)
                .hotspot(false)
                .start(createStartLeg(startChromosome, startPosition, type).build())
                .end(createEndLeg(endChromosome, endPosition, type).build())
                .startContext(dummyContext());
    }

    public static ImmutableStructuralVariantImpl.Builder createStructuralVariantSingleBreakend(
            final String startChromosome, final int startPosition, double startVaf)
    {
        return ImmutableStructuralVariantImpl.builder()
                .id(Strings.EMPTY)
                .insertSequence(Strings.EMPTY)
                .insertSequenceAlignments(Strings.EMPTY)
                .qualityScore(0)
                .hotspot(false)
                .type(StructuralVariantType.BND)
                .start(createStartLeg(startChromosome, startPosition, StructuralVariantType.BND).alleleFrequency(startVaf).build())
                .startContext(dummyContext());
    }

    @NotNull
    public static ImmutableStructuralVariantLegImpl.Builder createStartLeg(final String startChromosome, final int startPosition,
            final StructuralVariantType type)
    {
        final byte startOrientation = switch(type)
        {
            case DUP -> -1;
            case BND, INV -> 1;
            default -> 1;
        };

        return ImmutableStructuralVariantLegImpl.builder()
                .chromosome(startChromosome)
                .position(startPosition)
                .homology(Strings.EMPTY)
                .orientation(startOrientation)
                .anchoringSupportDistance(0);
    }

    private static ImmutableStructuralVariantLegImpl.Builder createEndLeg(
            final String endChromosome, final int endPosition, final StructuralVariantType type)
    {
        final byte endOrientation = switch(type)
        {
            case DUP -> 1;
            case BND, INV -> 1;
            default -> -1;
        };

        return ImmutableStructuralVariantLegImpl.builder()
                .chromosome(endChromosome)
                .position(endPosition)
                .orientation(endOrientation)
                .homology("")
                .anchoringSupportDistance(0);
    }

    private static VariantContext dummyContext()
    {
        final Collection<Allele> alleles = Lists.newArrayList(Allele.create("N", true));

        return new VariantContextBuilder("purple", "2", 1, 1, alleles)
                .noGenotypes()
                .make();
    }
}
