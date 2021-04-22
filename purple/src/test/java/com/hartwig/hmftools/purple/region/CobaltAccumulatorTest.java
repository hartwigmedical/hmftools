package com.hartwig.hmftools.purple.region;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelector;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.purple.PurpleTestUtils;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CobaltAccumulatorTest {

    private static final String CHROMOSOME = "1";
    private static final int WINDOW_SIZE = 1000;

    @Test
    public void testIncludeInside() {
        assertTumorCount(2, 5001, 10000, 4001, 5001, 9001, 10001);
    }

    @Test
    public void testExcludePartialBefore() {
        assertTumorCount(1, 5301, 10000, 4001, 5001, 9001, 10001);
    }

    @Test
    public void testExcludePartialAfter() {
        assertTumorCount(1, 5001, 9999, 4001, 5001, 9001, 10001);
    }

    private void assertTumorCount(int expectedCount, long regionStart, long regionEnd, long... ratios) {
        final GenomeRegion region = GenomeRegions.create(CHROMOSOME, regionStart, regionEnd);
        final ObservedRegionFactory.CobaltAccumulator accumulator = new ObservedRegionFactory.CobaltAccumulator(WINDOW_SIZE, region);
        final GenomePositionSelector<CobaltRatio> selector = createSelector(ratios);
        selector.select(region, accumulator);
        assertEquals(expectedCount, accumulator.tumorCount());
    }

    @NotNull
    private static GenomePositionSelector<CobaltRatio> createSelector(long... ratioPositions) {
        List<CobaltRatio> ratios = Lists.newArrayList();
        for (long ratioPosition : ratioPositions) {
            ratios.add(ratio(ratioPosition));
        }

        return GenomePositionSelectorFactory.create(ratios);
    }

    @NotNull
    private static CobaltRatio ratio(long position) {
        return PurpleTestUtils.cobalt(CHROMOSOME, position, 1).build();
    }
}
