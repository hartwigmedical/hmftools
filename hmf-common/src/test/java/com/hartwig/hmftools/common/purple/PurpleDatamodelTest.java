package com.hartwig.hmftools.common.purple;

import com.hartwig.hmftools.common.purple.copynumber.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.region.ImmutableEnrichedRegion;
import com.hartwig.hmftools.common.purple.region.ImmutableFittedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegionStatus;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PurpleDatamodelTest {

    private static final String CHROMOSOME = "1";

    @Test
    public void testDefaultFittedRegion() {
        createDefaultFittedRegion(CHROMOSOME, 1, 100).build();
    }

    @Test
    public void testDefaultCopyNumber() {
        createCopyNumber(CHROMOSOME, 1, 100, 2).build();
    }

    @NotNull
    public static ImmutablePurpleCopyNumber.Builder createCopyNumber(@NotNull final String chromosome, final long start, final long end,
            final double copyNumber) {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome(chromosome)
                .start(start)
                .end(end)
                .averageTumorCopyNumber(copyNumber)
                .support(SegmentSupport.NONE)
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
                .observedTumorRatioCount(1)
                .observedNormalRatio(1)
                .ratioSupport(true)
                .status(ObservedRegionStatus.SOMATIC)
                .support(SegmentSupport.NONE);
    }

    @NotNull
    public static ImmutableFittedRegion.Builder createDefaultFittedRegion(@NotNull final String chromosome, final long start,
            final long end) {
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
                .modelPloidy(0)
                .modelBAF(0)
                .modelTumorRatio(0)
                .refNormalisedCopyNumber(2)
                .ratioSupport(true)
                .support(SegmentSupport.NONE)
                .bafDeviation(0);
    }


}
