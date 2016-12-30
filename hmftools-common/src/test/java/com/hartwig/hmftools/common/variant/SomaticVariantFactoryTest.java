package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class SomaticVariantFactoryTest {

    @Test
    public void incorrectAFFieldYieldsNaN() {
        final String line = "0 \t 1 \t 2 \t 3 \t 4 \t 5 \t 6 \t ;set=Intersection; \t 8 \t 0/1:60:113";

        final SomaticVariant variant = SomaticVariantFactory.fromVCFLine(line);
        assertEquals(Double.NaN, variant.alleleFrequency(), 1.0e-10);
    }

    @Test
    public void correctDBSNPAndCOSMIC() {
        final String both = "0 \t 1 \t rs:1;COSM2 \t 3 \t 4 \t 5 \t 6 \t ;set=Intersection; \t 8 \t 0/1:60:113";

        final SomaticVariant hasBoth = SomaticVariantFactory.fromVCFLine(both);
        assertTrue(hasBoth.isDBSNP());
        assertTrue(hasBoth.isCOSMIC());

        final String dbsnpOnly = "0 \t 1 \t rs:1 \t 3 \t 4 \t 5 \t 6 \t ;set=Intersection; \t 8 \t 0/1:60:113";
        final SomaticVariant hasDBSNPOnly = SomaticVariantFactory.fromVCFLine(dbsnpOnly);
        assertTrue(hasDBSNPOnly.isDBSNP());
        assertFalse(hasDBSNPOnly.isCOSMIC());

        final String cosmicOnly = "0 \t 1 \t COSM2 \t 3 \t 4 \t 5 \t 6 \t ;set=Intersection; \t 8 \t 0/1:60:113";
        final SomaticVariant hasCOSMICOnly = SomaticVariantFactory.fromVCFLine(cosmicOnly);
        assertFalse(hasCOSMICOnly.isDBSNP());
        assertTrue(hasCOSMICOnly.isCOSMIC());

        final String none = "0 \t 1 \t 2 \t 3 \t 4 \t 5 \t 6 \t ;set=Intersection; \t 8 \t 0/1:60:113";
        final SomaticVariant hasNone = SomaticVariantFactory.fromVCFLine(none);
        assertFalse(hasNone.isDBSNP());
        assertFalse(hasNone.isCOSMIC());
    }
}