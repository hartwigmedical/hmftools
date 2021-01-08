package com.hartwig.hmftools.common.genome.region;

import static com.google.common.collect.Lists.newArrayList;
import static com.hartwig.hmftools.common.genome.region.GenomeRegions.create;
import static com.hartwig.hmftools.common.genome.region.GenomeRegionsValidation.isSubset;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Test;

public class GenomeRegionsValidationTest {

    @Test
    public void testCombinedCoverage() {
        final List<GenomeRegion> superset = newArrayList(
                create("1", 100, 110),
                create("1", 111, 120),
                create("1", 122, 130),
                create("1", 131, 140)
        );
        assertFalse(isSubset(superset, newArrayList(create("1", 99, 101))));
        assertTrue(isSubset(superset, newArrayList(create("1", 100, 101))));
        assertTrue(isSubset(superset, newArrayList(create("1", 101, 101))));

        assertFalse(isSubset(superset, newArrayList(create("2", 101, 101))));


        assertFalse(isSubset(superset, newArrayList(create("1", 99, 120))));
        assertTrue(isSubset(superset, newArrayList(create("1", 100, 120))));
        assertFalse(isSubset(superset, newArrayList(create("1", 100, 121))));

        assertTrue(isSubset(superset, newArrayList(create("1", 122, 140))));
        assertFalse(isSubset(superset, newArrayList(create("1", 100, 140))));

    }

}
