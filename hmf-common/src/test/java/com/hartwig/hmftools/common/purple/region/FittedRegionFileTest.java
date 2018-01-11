package com.hartwig.hmftools.common.purple.region;

import static org.junit.Assert.assertEquals;

import java.util.Random;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FittedRegionFileTest {

    @Test
    public void testSupportOlderFileFormat() {
        final FittedRegion complete = createRandom(new Random());
        final FittedRegion expectedTruncated = ImmutableFittedRegion.builder().from(complete).ploidyPenalty(0).build();

        final String completeString = FittedRegionFile.toString(complete);
        final String values[] = completeString.split("\t");

        final StringJoiner truncatedStringJoiner = new StringJoiner("\t");
        for (int i = 0; i < values.length - 1; i++) {
            truncatedStringJoiner.add(values[i]);
        }

        final FittedRegion truncated = FittedRegionFile.fromString(truncatedStringJoiner.toString());
        assertEquals(expectedTruncated, truncated);
    }

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
                .modelPloidy(random.nextInt())
                .bafCount(random.nextInt())
                .observedBAF(random.nextDouble())
                .tumorBAF(random.nextDouble())
                .modelBAF(random.nextDouble())
                .bafDeviation(random.nextDouble())
                .observedTumorRatio(random.nextDouble())
                .observedNormalRatio(random.nextDouble())
                .modelTumorRatio(random.nextDouble())
                .cnvDeviation(random.nextDouble())
                .deviation(random.nextDouble())
                .tumorCopyNumber(random.nextDouble())
                .segmentTumorCopyNumber(random.nextDouble())
                .segmentBAF(random.nextDouble())
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
