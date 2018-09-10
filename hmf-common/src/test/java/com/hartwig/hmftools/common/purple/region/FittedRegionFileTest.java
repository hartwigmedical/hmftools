package com.hartwig.hmftools.common.purple.region;

import static org.junit.Assert.assertEquals;

import java.util.Random;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FittedRegionFileTest {

    @Test
    public void testToFromString() {
        final FittedRegion expected = createRandom(new Random());
        final FittedRegion decoded = FittedRegionFile.fromString(FittedRegionFile.toString(expected));
        assertEquals(expected, decoded);
    }

    @NotNull
    private static FittedRegion createRandom(@NotNull final Random random) {
        return ImmutableFittedRegion.builder()
                .chromosome("XYZ" + random.nextInt())
                .start(random.nextLong())
                .end(random.nextLong())
                .status(GermlineStatus.DIPLOID)
                .bafCount(random.nextInt())
                .observedBAF(nextDouble(random))
                .tumorBAF(nextDouble(random))
                .minorAllelePloidy(nextDouble(random))
                .minorAllelePloidyDeviation(nextDouble(random))
                .observedTumorRatio(nextDouble(random))
                .observedNormalRatio(nextDouble(random))
                .majorAllelePloidy(nextDouble(random))
                .majorAllelePloidyDeviation(nextDouble(random))
                .deviation(nextDouble(random))
                .tumorCopyNumber(nextDouble(random))
                .fittedTumorCopyNumber(nextDouble(random))
                .fittedBAF(nextDouble(random))
                .refNormalisedCopyNumber(nextDouble(random))
                .ratioSupport(random.nextBoolean())
                .support(SegmentSupport.BND)
                .depthWindowCount(random.nextInt())
                .svCluster(random.nextBoolean())
                .gcContent(nextDouble(random))
                .ploidyPenalty(nextDouble(random))
                .build();
    }

    private static double nextDouble(@NotNull final Random random) {
        return Math.round(random.nextDouble() * 10000) / 10000;
    }

}
