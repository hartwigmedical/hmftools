package com.hartwig.hmftools.common.purple;

import java.util.Random;

import com.hartwig.hmftools.common.cobalt.ImmutableCobaltRatio;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.copynumber.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.purity.ImmutableFittedPurity;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;
import com.hartwig.hmftools.common.purple.region.ImmutableEnrichedRegion;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;

public class PurpleTestUtils
{
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

}
