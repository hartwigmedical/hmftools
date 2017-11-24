package com.hartwig.hmftools.common.gene;

import static org.junit.Assert.assertEquals;

import java.util.Random;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GeneCopyNumberTest {

    @Test
    public void testCopyNumberZero() {
        final GeneCopyNumber copyNumber = createRandom(new Random());
        assertEquals(0, copyNumber.value());
    }

    @NotNull
    private static GeneCopyNumber createRandom(@NotNull Random random) {
        return ImmutableGeneCopyNumber.builder()
                .chromosome(String.valueOf(random.nextInt(22)))
                .start(random.nextLong())
                .end(random.nextLong())
                .gene("gene" + random.nextInt())
                .minCopyNumber(-2.3)
                .maxCopyNumber(random.nextDouble())
                .meanCopyNumber(random.nextDouble())
                .somaticRegions(random.nextInt())
                .germlineHomRegions(random.nextInt())
                .germlineHet2HomRegions(random.nextInt())
                .transcriptID("transcriptId" + random.nextInt())
                .transcriptVersion(random.nextInt())
                .chromosomeBand("chromosomeband" + random.nextInt())
                .build();
    }
}
