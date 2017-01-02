package com.hartwig.hmftools.patientreporter.slicing;

import static org.junit.Assert.assertEquals;

import java.util.SortedSet;

import com.google.common.collect.Sets;

import org.junit.Test;

public class GenomeRegionTest {

    @Test
    public void sortsAscending() {
        GenomeRegion first = new GenomeRegion(100, 200);
        GenomeRegion last = new GenomeRegion(300, 400);

        SortedSet<GenomeRegion> sorted = Sets.newTreeSet();
        sorted.add(last);
        sorted.add(first);

        assertEquals(first, sorted.first());
        assertEquals(last, sorted.last());
    }
}