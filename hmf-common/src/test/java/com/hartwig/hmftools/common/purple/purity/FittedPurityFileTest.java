package com.hartwig.hmftools.common.purple.purity;

import static com.hartwig.hmftools.common.purple.purity.FittedPurityFile.HEADER_PREFIX;
import static com.hartwig.hmftools.common.purple.purity.FittedPurityFile.fromLine;
import static com.hartwig.hmftools.common.purple.purity.FittedPurityFile.toLines;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Random;

import com.hartwig.hmftools.common.purple.gender.Gender;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FittedPurityFileTest {

    @Test
    public void testInputAndOutput() {
        final Random random = new Random();
        final PurityContext input = ImmutablePurityContext.builder()
                .version(random.nextInt() + "a")
                .bestFit(createRandomPurity(random))
                .score(createRandomScore(random))
                .gender(Gender.values()[random.nextInt(Gender.values().length)])
                .status(FittedPurityStatus.values()[random.nextInt(FittedPurityStatus.values().length)])
                .polyClonalProportion(random.nextDouble())
                .build();

        final List<String> lines = toLines(input);
        assertEquals(2, lines.size());
        assertTrue(lines.get(0).startsWith(HEADER_PREFIX));

        final PurityContext output = fromLine(lines.get(1));
        assertEquals(input, output);
    }

    @NotNull
    private static FittedPurityScore createRandomScore(@NotNull Random random) {
        return ImmutableFittedPurityScore.builder()
                .minPurity(random.nextDouble())
                .maxPurity(random.nextDouble())
                .minPloidy(random.nextDouble())
                .maxPloidy(random.nextDouble())
                .minDiploidProportion(random.nextDouble())
                .maxDiploidProportion(random.nextDouble())
                .build();
    }

    @NotNull
    static FittedPurity createRandomPurity(@NotNull Random random) {
        return ImmutableFittedPurity.builder()
                .purity(random.nextDouble())
                .normFactor(random.nextDouble())
                .score(random.nextDouble())
                .diploidProportion(random.nextDouble())
                .ploidy(random.nextDouble())
                .somaticDeviation(random.nextDouble())
                .build();
    }
}
