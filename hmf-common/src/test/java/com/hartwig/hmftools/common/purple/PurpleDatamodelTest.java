package com.hartwig.hmftools.common.purple;

import com.hartwig.hmftools.common.copynumber.freec.FreecStatus;
import com.hartwig.hmftools.common.copynumber.freec.ImmutableFreecGCContent;
import com.hartwig.hmftools.common.purple.copynumber.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.region.ImmutableEnrichedRegion;
import com.hartwig.hmftools.common.purple.region.ImmutableFittedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.purple.segment.StructuralVariantSupport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PurpleDatamodelTest {

    public static final String CHROMOSOME = "1";

    @Test
    public void testDefaultFittedRegion() {
        createDefaultFittedRegion(CHROMOSOME, 1, 100).build();
    }

    @Test
    public void testDefaultCopyNumber() {
        createCopyNumber(CHROMOSOME, 1, 100, 2).build();
    }

    @NotNull
    public static ImmutableFreecGCContent.Builder createGCContent(String chromosome, long position, double gcContent, double nonNPercentage,
            double mappablePercentage) {
        return ImmutableFreecGCContent.builder()
                .chromosome(chromosome)
                .position(position)
                .gcContent(gcContent)
                .nonNPercentage(nonNPercentage)
                .mappablePercentage(mappablePercentage);
    }

    @NotNull
    public static ImmutablePurpleCopyNumber.Builder createCopyNumber(@NotNull String chromosome, long start, long end, double copyNumber) {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome(chromosome)
                .start(start)
                .end(end)
                .averageTumorCopyNumber(copyNumber)
                .structuralVariantSupport(StructuralVariantSupport.NONE)
                .ratioSupport(true)
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
                .observedNormalRatio(1)
                .observedGCContent(1)
                .observedNonNPercentage(1)
                .observedMappablePercentage(1)
                .ratioSupport(true)
                .structuralVariantSupport(StructuralVariantSupport.NONE);
    }

    @NotNull
    public static ImmutableFittedRegion.Builder createDefaultFittedRegion(@NotNull final String chromosome, long start, long end) {
        final ObservedRegion observedRegion = createObservedRegion(chromosome, start, end).build();
        return ImmutableFittedRegion.builder()
                .from(observedRegion)
                .tumorCopyNumber(2)
                .broadBAF(0)
                .broadTumorCopyNumber(0)
                .segmentBAF(0)
                .segmentTumorCopyNumber(0)
                .cnvDeviation(0)
                .deviation(0)
                .fittedPloidy(0)
                .modelBAF(0)
                .modelTumorRatio(0)
                .status(FreecStatus.UNKNOWN)
                .refNormalisedCopyNumber(2)
                .ratioSupport(true)
                .structuralVariantSupport(StructuralVariantSupport.NONE)
                .bafDeviation(0);
    }

}
