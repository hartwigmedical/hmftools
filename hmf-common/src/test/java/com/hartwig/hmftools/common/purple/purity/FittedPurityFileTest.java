package com.hartwig.hmftools.common.purple.purity;

import static com.hartwig.hmftools.common.purple.purity.FittedPurityFile.DELIMITER;
import static com.hartwig.hmftools.common.purple.purity.FittedPurityFile.HEADER_PREFIX;
import static com.hartwig.hmftools.common.purple.purity.FittedPurityFile.bestFit;
import static com.hartwig.hmftools.common.purple.purity.FittedPurityFile.gender;
import static com.hartwig.hmftools.common.purple.purity.FittedPurityFile.polyClonalProportion;
import static com.hartwig.hmftools.common.purple.purity.FittedPurityFile.score;
import static com.hartwig.hmftools.common.purple.purity.FittedPurityFile.status;
import static com.hartwig.hmftools.common.purple.purity.FittedPurityFile.toLines;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.gender.Gender;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FittedPurityFileTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void testInputAndOutput() {
        final Random random = new Random();
        final PurityContext context = ImmutablePurityContext.builder()
                .bestFit(createRandomPurity(random))
                .bestPerPurity(Lists.newArrayList())
                .score(createRandomScore(random))
                .gender(Gender.values()[random.nextInt(Gender.values().length)])
                .status(FittedPurityStatus.values()[random.nextInt(FittedPurityStatus.values().length)])
                .polyClonalProportion(random.nextDouble())
                .build();

        final List<String> lines = toLines(context);
        assertEquals(2, lines.size());
        assertTrue(lines.get(0).startsWith(HEADER_PREFIX));

        final String[] values = lines.get(1).split(DELIMITER);
        assertEquals(context.bestFit(), bestFit(values));
        assertEquals(context.score(), score(values));
        assertEquals(context.gender(), gender(values));
        assertEquals(context.status(), status(values));
        assertEquals(context.polyClonalProportion(), polyClonalProportion(values), EPSILON);
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
                .build();
    }
}
