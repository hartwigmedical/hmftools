package com.hartwig.hmftools.common.purple.purity;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FittedPurityFileTest {

    @Test
    public void testHeaderIsGenerated() {
        int size = 4;
        final List<FittedPurity> purities = createPurities(size);
        final List<String> toLines = FittedPurityFile.toLines(purities);
        assertEquals(size + 1, toLines.size());
        assertEquals('#', toLines.get(0).charAt(0));
    }

    @Test
    public void testInputAndOutput() {
        final List<FittedPurity> expected = createPurities(5);
        final List<FittedPurity> victim = FittedPurityFile.fromLines(FittedPurityFile.toLines(expected));

        assertEquals(expected.size(), victim.size());
        for (int i = 0; i < expected.size(); i++) {
            assertEquals(expected.get(i), victim.get(i));
        }
    }

    @NotNull
    private List<FittedPurity> createPurities(int count) {
        Random random = new Random();
        final List<FittedPurity> result = Lists.newArrayList();
        for (int i = 0; i < count; i++) {
            result.add(createRandom(random));
        }
        return result;
    }

    @NotNull
    private FittedPurity createRandom(@NotNull Random random) {
        return ImmutableFittedPurity.builder()
                .purity(random.nextDouble())
                .normFactor(random.nextDouble())
                .score(random.nextDouble())
                .modelBAFDeviation(random.nextDouble())
                .diploidProportion(random.nextDouble())
                .build();
    }
}
