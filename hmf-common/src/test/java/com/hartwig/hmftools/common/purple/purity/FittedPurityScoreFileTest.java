package com.hartwig.hmftools.common.purple.purity;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Random;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FittedPurityScoreFileTest {

    @Test
    public void testHeaderIsGenerated() {
        final FittedPurityScore score = createRandom(new Random());
        final List<String> toLines = FittedPurityScoreFile.toLines(score);
        assertEquals(2, toLines.size());
        assertTrue(toLines.get(0).startsWith(FittedPurityScoreFile.HEADER_PREFIX));
    }

    @Test
    public void testInputAndOutput() {
        final FittedPurityScore expected = createRandom(new Random());
        final FittedPurityScore victim = FittedPurityScoreFile.fromLines(FittedPurityScoreFile.toLines(expected));

        assertEquals(expected, victim);
    }

    @NotNull
    private FittedPurityScore createRandom(@NotNull Random random) {
        return ImmutableFittedPurityScore.builder()
                .polyclonalProportion(random.nextDouble())
                .minPurity(random.nextDouble())
                .maxPurity(random.nextDouble())
                .minPloidy(random.nextDouble())
                .maxPloidy(random.nextDouble())
                .build();
    }
}
