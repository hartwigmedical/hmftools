package com.hartwig.hmftools.common.purple.copynumber;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class ExtendDiploidBAFTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void testNeighbourIsHetrozygousDiploid() {
        assertBAF(1.0, 0.9, 0.5, 2);
        assertBAF(1.0, 1.0, 0.5, 2);
        assertBAF(2/3d, 1.5, 0.5, 2);
        assertBAF(0.5, 2.0, 0.5, 2);
        assertBAF(2/3d, 3.0, 0.5, 2);
        assertBAF(0.75, 4.0, 0.5, 2);
        assertBAF(0.80, 5.0, 0.5, 2);
    }

    @Test
    public void testNeighbourIsHomozygous() {
        assertBAF(1.0, 0.9, 1, 1);
        assertBAF(1.0, 1.0, 1, 1);
        assertBAF(1.0, 1.5, 1, 2);
        assertBAF(1.0, 2.0, 1, 2);
        assertBAF(1.0, 3.0, 1, 3);
        assertBAF(1.0, 4.0, 1, 4);
        assertBAF(1.0, 5.0, 1, 5);
    }

    @Test
    public void testNeighbourHasCN6() {
        assertBAF(1d, 2, 0.5, 6);
        assertBAF(1d, 3, 0.5, 6);
        assertBAF(3/4d, 4, 0.5, 6);
        assertBAF(3/5d, 5, 0.5, 6);

        assertBAF(1d, 2, 2/3d, 6);
        assertBAF(2/3d, 3, 2/3d, 6);
        assertBAF(2/4d, 4, 2/3d, 6);
        assertBAF(3/5d, 5, 2/3d, 6);
    }

    private void assertBAF(double expectedBAF, double regionCopyNumber, double neighbourBAF, double neighbourCopyNumber) {
        assertEquals(expectedBAF, ExtendDiploidBAF.estimateBAF(regionCopyNumber, neighbourBAF, neighbourCopyNumber), EPSILON);
    }

}
