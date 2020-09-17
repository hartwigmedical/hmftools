package com.hartwig.hmftools.common.genome.chromosome;

import com.hartwig.hmftools.common.cobalt.ImmutableMedianRatio;
import com.hartwig.hmftools.common.cobalt.MedianRatio;

import org.junit.Test;

public class CobaltChromosomesTest {

    @Test
    public void testMosiacX() {

    }


    private static MedianRatio create(String contig, double ratio) {
        return create(contig, ratio, CobaltChromosomes.MIN_RATIO_COUNT);
    }

    private static MedianRatio create(String contig, double ratio, int count) {
        return ImmutableMedianRatio.builder().count(count).chromosome(contig).medianRatio(ratio).build();
    }

}
