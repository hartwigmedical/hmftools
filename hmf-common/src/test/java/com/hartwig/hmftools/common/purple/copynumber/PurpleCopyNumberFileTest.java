package com.hartwig.hmftools.common.purple.copynumber;

import static com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile.HEADER_PREFIX;
import static com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile.fromLines;
import static com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile.toLines;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PurpleCopyNumberFileTest {

    @Test
    public void testHeaderIsGenerated() {
        int size = 4;
        final List<PurpleCopyNumber> copyNumbers = create(size);
        final List<String> toLines = toLines(copyNumbers);
        assertEquals(size + 1, toLines.size());
        assertTrue(toLines.get(0).startsWith(HEADER_PREFIX));
    }

    @Test
    public void testInputAndOutput() {
        final List<PurpleCopyNumber> expected = create(5);
        final List<PurpleCopyNumber> victim = fromLines(toLines(expected));

        assertEquals(expected.size(), victim.size());
        for (int i = 0; i < expected.size(); i++) {
            assertEquals(expected.get(i), victim.get(i));
        }
    }

    @NotNull
    private List<PurpleCopyNumber> create(int count) {
        Random random = new Random();
        final List<PurpleCopyNumber> result = Lists.newArrayList();
        for (int i = 0; i < count; i++) {
            result.add(createRandom(random));
        }
        return result;
    }

    @NotNull
    private PurpleCopyNumber createRandom(@NotNull Random random) {
        return PurpleDatamodelTest.createCopyNumber(random.nextInt(22) + "", random.nextLong(), random.nextLong(), random.nextDouble())
                .bafCount(random.nextInt())
                .averageObservedBAF(random.nextDouble())
                .averageActualBAF(random.nextDouble())
                .segmentStartSupport(SegmentSupport.values()[random.nextInt(SegmentSupport.values().length)])
                .segmentEndSupport(SegmentSupport.values()[random.nextInt(SegmentSupport.values().length)])
                .method(CopyNumberMethod.values()[random.nextInt(CopyNumberMethod.values().length)])
                .build();
    }
}
