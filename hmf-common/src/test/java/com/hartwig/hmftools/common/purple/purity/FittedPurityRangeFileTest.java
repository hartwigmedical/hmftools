package com.hartwig.hmftools.common.purple.purity;

import static com.hartwig.hmftools.common.purple.PurpleTestUtils.createRandomPurityBuilder;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FittedPurityRangeFileTest {

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

    @NotNull
    private static List<FittedPurity> create(int count) {
        Random random = new Random();
        final List<FittedPurity> result = Lists.newArrayList();
        for (int i = 0; i < count; i++) {
            result.add(createRandomPurityBuilder(random).purity(i).build());
        }
        return result;
    }
}
