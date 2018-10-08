package com.hartwig.hmftools.common.variant.recovery;

import static org.junit.Assert.assertEquals;

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

}
