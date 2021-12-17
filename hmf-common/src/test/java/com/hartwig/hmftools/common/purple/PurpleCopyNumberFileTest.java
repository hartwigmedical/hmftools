package com.hartwig.hmftools.common.purple;

import static com.hartwig.hmftools.common.purple.PurpleTestUtils.createCopyNumber;
import static com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile.fromLines;
import static com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile.toLines;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.nio.charset.Charset;
import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
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
        assertTrue(toLines.get(0).startsWith("chr"));
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
    private static List<PurpleCopyNumber> create(int count) {
        Random random = new Random();
        final List<PurpleCopyNumber> result = Lists.newArrayList();
        for (int i = 0; i < count; i++) {
            result.add(createRandom(random));
        }
        return result;
    }

    @NotNull
    private static PurpleCopyNumber createRandom(@NotNull Random random) {
        return createCopyNumber(random.nextInt(22) + "", random.nextInt(), random.nextInt(), nextDouble(random))
                .bafCount(random.nextInt())
                .averageObservedBAF(nextDouble(random))
                .averageActualBAF(nextDouble(random))
                .segmentStartSupport(SegmentSupport.values()[random.nextInt(SegmentSupport.values().length)])
                .segmentEndSupport(SegmentSupport.values()[random.nextInt(SegmentSupport.values().length)])
                .method(CopyNumberMethod.values()[random.nextInt(CopyNumberMethod.values().length)])
                .depthWindowCount(random.nextInt())
                .gcContent(nextDouble(random))
                .minStart(random.nextLong())
                .maxStart(random.nextLong())
                .build();
    }

    private static double nextDouble(@NotNull final Random random) {
        return Math.round(random.nextDouble() * 10000D) / 10000D;
    }
}
