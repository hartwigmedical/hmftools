package com.hartwig.hmftools.common.purple.purity;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FittedPurityFileTest {

    @Test
    public void testHeaderIsGenerated() {
        int size = 4;
        final List<FittedPurity> purities = create(size);
        final List<String> toLines = FittedPurityFile.toLines(purities);
        assertEquals(size + 1, toLines.size());
        assertTrue(toLines.get(0).startsWith(FittedPurityFile.HEADER_PREFIX));
    }

    @Test
    public void testInputAndOutput() {
        final List<FittedPurity> expected = create(5);
        final List<FittedPurity> victim = FittedPurityFile.fromLines(FittedPurityFile.toLines(expected));

        assertEquals(expected.size(), victim.size());
        for (int i = 0; i < expected.size(); i++) {
            assertEquals(expected.get(i), victim.get(i));
        }
    }

    @NotNull
    private List<FittedPurity> create(int count) {
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
