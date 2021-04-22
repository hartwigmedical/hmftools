package com.hartwig.hmftools.purple.purity;

import static com.hartwig.hmftools.common.purple.PurpleTestUtils.createRandomPurityBuilder;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.purity.BestFit;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class BestFitTest {

    private final Random random = new Random();

    @Test
    public void testBestFitPerPurity() {
        final FittedPurity fp1 = createRandomPurity(0.3, 0.3, 1.9);
        final FittedPurity fp2 = createRandomPurity(0.3, 0.2, 2.3);
        final FittedPurity fp3 = createRandomPurity(0.4, 0.4, 2);
        final FittedPurity fp4 = createRandomPurity(0.4, 0.3, 2);

        final List<FittedPurity> all = Lists.newArrayList(fp1, fp2, fp3, fp4);
        Collections.shuffle(all);

        final List<FittedPurity> bestFitPerPurity = BestFit.bestFitPerPurity(all);
        assertEquals(2, bestFitPerPurity.size());
        assertEquals(fp2, bestFitPerPurity.get(0));
        assertEquals(fp4, bestFitPerPurity.get(1));
    }

    @Test
    public void testMostDiploidPurity() {
        final FittedPurity fp1 = createRandomPurity(0.3, 0.3, 2.3);
        final FittedPurity fp2 = createRandomPurity(0.3, 0.2, 1.9);
        final FittedPurity fp3 = createRandomPurity(0.4, 0.4, 1.8);
        final FittedPurity fp4 = createRandomPurity(0.4, 0.3, 2.05);

        final List<FittedPurity> all = Lists.newArrayList(fp1, fp2, fp3, fp4);
        Collections.shuffle(all);

        final List<FittedPurity> result = BestFit.mostDiploidPerPurity(all);

        assertEquals(2, result.size());
        assertEquals(fp2, result.get(0));
        assertEquals(fp4, result.get(1));
    }

    @NotNull
    private FittedPurity createRandomPurity(double purity, double score, double ploidy)
    {
        return createRandomPurityBuilder(random).purity(purity).score(score).ploidy(ploidy).build();
    }

}
