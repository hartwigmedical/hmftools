package com.hartwig.hmftools.common.region;

import static org.junit.Assert.assertEquals;

import java.util.SortedSet;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.region.bed.ImmutableBEDGenomeRegion;

import org.junit.Test;

public class GenomeRegionTest {

    @Test
    public void sortsAscending() {
        final GenomeRegion first = ImmutableBEDGenomeRegion.of("X", 100, 200);
        final GenomeRegion last = ImmutableBEDGenomeRegion.of("X", 300, 400);

        final SortedSet<GenomeRegion> sorted = Sets.newTreeSet();
        sorted.add(last);
        sorted.add(first);

        assertEquals(first, sorted.first());
        assertEquals(last, sorted.last());
    }
}