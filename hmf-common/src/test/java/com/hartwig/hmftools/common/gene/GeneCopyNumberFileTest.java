package com.hartwig.hmftools.common.gene;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GeneCopyNumberFileTest {


    @Test
    public void testInputAndOutput() {
        final List<GeneCopyNumber> expected = create(5);
        final List<GeneCopyNumber> victim = GeneCopyNumberFile.fromLines(GeneCopyNumberFile.toLines(expected));

        assertEquals(expected.size(), victim.size());
        for (int i = 0; i < expected.size(); i++) {
            assertEquals(expected.get(i), victim.get(i));
        }
    }

    @NotNull
    private List<GeneCopyNumber> create(int count) {
        Random random = new Random();
        final List<GeneCopyNumber> result = Lists.newArrayList();
        for (int i = 0; i < count; i++) {
            result.add(createRandom(random));
        }
        return result;
    }

    @NotNull
    private static GeneCopyNumber createRandom(@NotNull Random random) {
        return ImmutableGeneCopyNumber.builder()
                .chromosome(String.valueOf(random.nextInt(22)))
                .start(random.nextLong())
                .end(random.nextLong())
                .gene("gene" + random.nextInt())
                .minCopyNumber(random.nextDouble())
                .maxCopyNumber(random.nextDouble())
                .meanCopyNumber(random.nextDouble())
                .regions(random.nextInt())
                .transcriptID("transcriptId" + random.nextInt())
                .transcriptVersion(random.nextInt())
                .chromosomeBand("chromosomeband" + random.nextInt())
                .build();
    }
}
