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
                .observedBAF(random.nextDouble())
                .tumorBAF(random.nextDouble())
                .minorAllelePloidy(random.nextDouble())
                .minorAllelePloidyDeviation(random.nextDouble())
                .observedTumorRatio(random.nextDouble())
                .observedNormalRatio(random.nextDouble())
                .majorAllelePloidy(random.nextDouble())
                .majorAllelePloidyDeviation(random.nextDouble())
                .deviation(random.nextDouble())
                .tumorCopyNumber(random.nextDouble())
                .fittedTumorCopyNumber(random.nextDouble())
                .fittedBAF(random.nextDouble())
                .refNormalisedCopyNumber(random.nextDouble())
                .ratioSupport(random.nextBoolean())
                .support(SegmentSupport.BND)
                .observedTumorRatioCount(random.nextInt())
                .svCluster(random.nextBoolean())
                .gcContent(random.nextDouble())
                .ploidyPenalty(random.nextDouble())
                .build();
    }
}
