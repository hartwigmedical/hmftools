package com.hartwig.hmftools.common.gene;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GeneCopyNumberFileTest {

    @Test
    public void testSupportOlderFileFormat() {
        final String oldVersion =
                "1\t11869\t14409\tDDX11L1\t2.60025049329193\t2.60025049329193\t2.60025049329193\t1\t0\t0\tENST00000456328\t2\tp36.33";

        GeneCopyNumber copyNumber = GeneCopyNumberFile.fromString(oldVersion);
        assertEquals(11869, copyNumber.start());
        assertEquals("p36.33", copyNumber.chromosomeBand());
    }

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
            result.add(createRandom(random).build());
        }
        return result;
    }

    @NotNull
    static ImmutableGeneCopyNumber.Builder createRandom(@NotNull Random random) {
        return ImmutableGeneCopyNumber.builder()
                .chromosome(String.valueOf(random.nextInt(22)))
                .start(random.nextLong())
                .end(random.nextLong())
                .gene("gene" + random.nextInt())
                .minCopyNumber(nextDouble(random))
                .maxCopyNumber(nextDouble(random))
                .somaticRegions(random.nextInt())
                .germlineHomRegions(random.nextInt())
                .germlineHet2HomRegions(random.nextInt())
                .transcriptID("transcriptId" + random.nextInt())
                .transcriptVersion(random.nextInt())
                .chromosomeBand("chromosomeband" + random.nextInt())
                .minRegions(random.nextInt())
                .minRegionStart(random.nextLong())
                .minRegionEnd(random.nextLong())
                .minRegionMethod(randomMethod(random))
                .minRegionStartSupport(randomSupport(random))
                .minRegionEndSupport(randomSupport(random))
                .nonsenseBiallelicCount(random.nextInt())
                .nonsenseNonBiallelicCount(random.nextInt())
                .nonsenseNonBiallelicPloidy(nextDouble(random))
                .spliceBiallelicCount(random.nextInt())
                .spliceNonBiallelicCount(random.nextInt())
                .spliceNonBiallelicPloidy(nextDouble(random))
                .missenseBiallelicCount(random.nextInt())
                .missenseNonBiallelicCount(random.nextInt())
                .missenseNonBiallelicPloidy(random.nextInt())
                .minMinorAllelePloidy(nextDouble(random));
    }

    @NotNull
    private static CopyNumberMethod randomMethod(@NotNull Random random) {
        return CopyNumberMethod.values()[random.nextInt(CopyNumberMethod.values().length)];
    }

    @NotNull
    private static SegmentSupport randomSupport(@NotNull Random random) {
        return SegmentSupport.values()[random.nextInt(SegmentSupport.values().length)];
    }
    
    private static double nextDouble(@NotNull final Random random) {
        return Math.round(random.nextDouble() * 10000) / 10000;
    }
    

}
