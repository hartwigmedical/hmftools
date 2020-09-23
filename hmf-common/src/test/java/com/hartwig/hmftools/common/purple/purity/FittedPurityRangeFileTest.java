package com.hartwig.hmftools.common.purple.purity;

import static com.hartwig.hmftools.common.purple.purity.PurityContextFileTest.createRandomPurityBuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

public class FittedPurityRangeFileTest {

    private Random random;

    @Before
    public void setup() {
        random = new Random();
    }

    @Test
    public void testHeaderIsGenerated() {
        int size = 4;
        final List<FittedPurity> purities = create(size);
        final List<String> toLines = FittedPurityRangeFile.toLines(purities);
        assertEquals(size + 1, toLines.size());
        assertTrue(toLines.get(0).startsWith("purity"));
    }

    @Test
    public void testInputAndOutput() {
        final List<FittedPurity> expected = create(5).stream().sorted().collect(Collectors.toList());
        final List<FittedPurity> victim = FittedPurityRangeFile.fromLines(FittedPurityRangeFile.toLines(expected));

        assertEquals(expected.size(), victim.size());
        for (int i = 0; i < expected.size(); i++) {
            assertEquals(expected.get(i), victim.get(i));
        }
    }

    @Test
    public void testBestFitPerPurity() {
        final FittedPurity fp1 = createRandomPurity(0.3, 0.3);
        final FittedPurity fp2 = createRandomPurity(0.3, 0.2);
        final FittedPurity fp3 = createRandomPurity(0.4, 0.4);
        final FittedPurity fp4 = createRandomPurity(0.4, 0.3);

        final List<FittedPurity> all = Lists.newArrayList(fp1, fp2, fp3, fp4);
        Collections.shuffle(all);

        final List<FittedPurity> bestFitPerPurity = FittedPurityRangeFile.bestFitPerPurity(all);
        assertEquals(2, bestFitPerPurity.size());
        assertEquals(fp2, bestFitPerPurity.get(0));
        assertEquals(fp4, bestFitPerPurity.get(1));
    }

    @NotNull
    private static List<FittedPurity> create(int count) {
        Random random = new Random();
        final List<FittedPurity> result = Lists.newArrayList();
        for (int i = 0; i < count; i++) {
            result.add(createRandomPurityBuilder(random).purity(i).build());
        }
        return result;
    }

    @NotNull
    private FittedPurity createRandomPurity(double purity, double score) {
        return createRandomPurityBuilder(random).purity(purity).score(score).build();
    }
}
