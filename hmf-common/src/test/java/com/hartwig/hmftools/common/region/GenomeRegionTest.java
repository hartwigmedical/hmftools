package com.hartwig.hmftools.common.region;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.region.bed.ImmutableBEDGenomeRegion;
import org.junit.Test;

import java.util.SortedSet;

import static org.junit.Assert.assertEquals;

public class GenomeRegionTest {

    @Test
    public void sortsAscending() {
        final GenomeRegion first = ImmutableBEDGenomeRegion.of("X", 100, 200, null);
        final GenomeRegion last = ImmutableBEDGenomeRegion.of("X", 300, 400, null);

        final SortedSet<GenomeRegion> sorted = Sets.newTreeSet();
        sorted.add(last);
        sorted.add(first);

        assertEquals(first, sorted.first());
        assertEquals(last, sorted.last());
    }
}