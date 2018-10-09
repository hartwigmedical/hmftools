package com.hartwig.hmftools.common.variant.recovery;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionFactory;

import org.junit.Test;

public class StructuralVariantRecoveryTest {

    @Test
    public void testOrientation() {
        assertEquals(1, StructuralVariantRecovery.orientation("C[17:59493156["));
        assertEquals(-1, StructuralVariantRecovery.orientation("]17:59493156]C"));
    }

    @Test
    public void testMate() {
        assertEquals("17:59493156", StructuralVariantRecovery.mateLocation("C[17:59493156["));
        assertEquals("17:59493156", StructuralVariantRecovery.mateLocation("]17:59493156]C"));
    }

    @Test
    public void testClosestRegion() {
        final GenomeRegion start = GenomeRegionFactory.create("1", 1, 10000);
        final GenomeRegion middle1 = GenomeRegionFactory.create("1", 10001, 20000);
        final GenomeRegion middle2 = GenomeRegionFactory.create("1", 20001, 30000);
        final GenomeRegion end = GenomeRegionFactory.create("1", 30001, 40000);
        final List<GenomeRegion> regions = Lists.newArrayList(start, middle1, middle2, end);

        assertEquals(start, StructuralVariantRecovery.closest(1, regions));
        assertEquals(start, StructuralVariantRecovery.closest(5001, regions));
        assertEquals(middle1, StructuralVariantRecovery.closest(5002, regions));
        assertEquals(middle1, StructuralVariantRecovery.closest(15001, regions));
        assertEquals(middle2, StructuralVariantRecovery.closest(15002, regions));
        assertEquals(middle2, StructuralVariantRecovery.closest(25001, regions));
        assertEquals(end, StructuralVariantRecovery.closest(25002, regions));
        assertEquals(end, StructuralVariantRecovery.closest(40000, regions));
    }

}
