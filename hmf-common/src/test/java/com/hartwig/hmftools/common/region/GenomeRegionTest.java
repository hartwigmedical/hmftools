package com.hartwig.hmftools.common.region;

import static org.junit.Assert.assertEquals;

import java.util.SortedSet;

import com.google.common.collect.Sets;

import org.junit.Test;

public class GenomeRegionTest {

    @Test
    public void sortsAscending() {
        final GenomeRegion first = GenomeRegionFactory.create("X", 100, 200);
        final GenomeRegion last = GenomeRegionFactory.create("X", 300, 400);

        final SortedSet<GenomeRegion> sorted = Sets.newTreeSet();
        sorted.add(last);
        sorted.add(first);

        assertEquals(first, sorted.first());
        assertEquals(last, sorted.last());
    }
}