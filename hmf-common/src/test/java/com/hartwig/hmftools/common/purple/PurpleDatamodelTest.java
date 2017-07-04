package com.hartwig.hmftools.common.purple;

import com.hartwig.hmftools.common.copynumber.freec.FreecStatus;
import com.hartwig.hmftools.common.purple.copynumber.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.region.ImmutableFittedRegion;
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
    public static ImmutablePurpleCopyNumber.Builder createCopyNumber(@NotNull String chromosome, long start, long end, double copyNumber) {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome(chromosome)
                .start(start)
                .end(end)
                .averageTumorCopyNumber(copyNumber)
                .bafCount(0)
                .averageObservedBAF(0.5)
                .averageActualBAF(0.5);
    }

    @NotNull
    public static ImmutableFittedRegion.Builder createDefaultFittedRegion(@NotNull String chromosome, long start, long end) {
        return ImmutableFittedRegion.builder()
                .chromosome(chromosome)
                .start(start)
                .end(end)
                .bafCount(1)
                .observedBAF(0.5)
                .tumorCopyNumber(2)
                .broadBAF(0)
                .broadTumorCopyNumber(0)
                .segmentBAF(0)
                .segmentTumorCopyNumber(0)
                .observedNormalRatio(1.0)
                .observedNormalRatio(1.0)
                .cnvDeviation(0)
                .deviation(0)
                .fittedPloidy(0)
                .modelBAF(0)
                .observedTumorRatio(0)
                .modelTumorRatio(0)
                .status(FreecStatus.UNKNOWN)
                .refNormalisedCopyNumber(2)
                .ratioSupport(true)
                .structuralVariantSupport(StructuralVariantSupport.NONE)
                .bafDeviation(0);
    }

}
