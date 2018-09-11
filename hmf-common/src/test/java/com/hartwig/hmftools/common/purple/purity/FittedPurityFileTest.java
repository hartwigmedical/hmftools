package com.hartwig.hmftools.common.purple.purity;

import static com.hartwig.hmftools.common.purple.purity.FittedPurityFile.HEADER_PREFIX;
import static com.hartwig.hmftools.common.purple.purity.FittedPurityFile.fromLine;
import static com.hartwig.hmftools.common.purple.purity.FittedPurityFile.toLines;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.nio.charset.Charset;
import java.util.List;
import java.util.Random;

import com.google.common.io.Resources;
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
                .polyClonalProportion(nextDouble(random))
                .build();

        final List<String> lines = toLines(input);
        assertEquals(2, lines.size());
        assertTrue(lines.get(0).startsWith(HEADER_PREFIX));

        final PurityContext output = fromLine(lines.get(1));
        assertEquals(input, output);
    }

    @Test
    public void testCompatibilityWith2_14() throws IOException {
        FittedPurityFile.fromLine(Resources.readLines(Resources.getResource("purple/v2-14.purple.purity"), Charset.defaultCharset())
                .get(1));
    }

    @NotNull
    private static FittedPurityScore createRandomScore(@NotNull Random random) {
        return ImmutableFittedPurityScore.builder()
                .minPurity(nextDouble(random))
                .maxPurity(nextDouble(random))
                .minPloidy(nextDouble(random))
                .maxPloidy(nextDouble(random))
                .minDiploidProportion(nextDouble(random))
                .maxDiploidProportion(nextDouble(random))
                .build();
    }

    @NotNull
    static FittedPurity createRandomPurity(@NotNull Random random) {
        return ImmutableFittedPurity.builder()
                .purity(nextDouble(random))
                .normFactor(nextDouble(random))
                .score(nextDouble(random))
                .diploidProportion(nextDouble(random))
                .ploidy(nextDouble(random))
                .somaticDeviation(nextDouble(random))
                .build();
    }

    private static double nextDouble(@NotNull final Random random) {
        return Math.round(random.nextDouble() * 10000) / 10000;
    }
}
