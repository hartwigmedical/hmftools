package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class SomaticVariantFactoryTest {

    private static final double EPSILON = 1.0e-10;

    @Test
    public void canReadCorrectSomaticVariant() {
        final String line = "15 \t 12345678 \t rs:1;UCSC \t C \t A,G \t <qual> \t <filter> "
                + " \t set=varscan-freebayes; \t <format> \t 0/1:60,60:121";

        final SomaticVariant variant = SomaticVariantFactory.fromVCFLine(line);
        assertEquals("15", variant.chromosome());
        assertEquals(12345678, variant.position());
        assertEquals(VariantType.SNP, variant.type());
        assertTrue(variant.isDBSNP());
        assertFalse(variant.isCOSMIC());

        assertEquals(2, variant.callerCount());
        assertTrue(variant.callers().contains(SomaticVariantConstants.FREEBAYES));
        assertTrue(variant.callers().contains(SomaticVariantConstants.VARSCAN));
        assertFalse(variant.callers().contains(SomaticVariantConstants.STRELKA));
        assertFalse(variant.callers().contains(SomaticVariantConstants.MUTECT));

        assertEquals(0.5, variant.alleleFrequency(), EPSILON);
    }

    @Test
    public void incorrectAFFieldYieldsNaN() {
        final String missingAFLine = "0 \t 1 \t 2 \t 3 \t 4 \t 5 \t 6 \t set=Intersection; \t 8 \t 9";
        final SomaticVariant missingAFVariant = SomaticVariantFactory.fromVCFLine(missingAFLine);
        assertEquals(Double.NaN, missingAFVariant.alleleFrequency(), EPSILON);

        final String missingRefCovLine = "0 \t 1 \t 2 \t 3 \t 4 \t 5 \t 6 \t set=Intersection; \t 8 \t 0/1:60:113";
        final SomaticVariant missingRefCovVariant = SomaticVariantFactory.fromVCFLine(missingRefCovLine);
        assertEquals(Double.NaN, missingRefCovVariant.alleleFrequency(), EPSILON);
    }

    @Test
    public void recognizeFilterInVarscan() {
        final String line = "0 \t 1 \t 2 \t 3 \t 4 \t 5 \t 6 \t set=freebayes-filterInVarscan; \t 8 \t 9";
        final SomaticVariant variant = SomaticVariantFactory.fromVCFLine(line);

        assertEquals(1, variant.callerCount());
        assertFalse(variant.callers().contains(SomaticVariantConstants.VARSCAN));
    }

    @Test
    public void correctDBSNPAndCOSMIC() {
        final String both = "0 \t 1 \t rs:1;COSM2 \t 3 \t 4 \t 5 \t 6 \t 7 \t 8 \t 9";

        final SomaticVariant hasBoth = SomaticVariantFactory.fromVCFLine(both);
        assertTrue(hasBoth.isDBSNP());
        assertTrue(hasBoth.isCOSMIC());

        final String dbsnpOnly = "0 \t 1 \t rs:1 \t 3 \t 4 \t 5 \t 6 \t 7 \t 8 \t 9";
        final SomaticVariant hasDBSNPOnly = SomaticVariantFactory.fromVCFLine(dbsnpOnly);
        assertTrue(hasDBSNPOnly.isDBSNP());
        assertFalse(hasDBSNPOnly.isCOSMIC());

        final String cosmicOnly = "0 \t 1 \t COSM2 \t 3 \t 4 \t 5 \t 6 \t 7 \t 8 \t 9";
        final SomaticVariant hasCOSMICOnly = SomaticVariantFactory.fromVCFLine(cosmicOnly);
        assertFalse(hasCOSMICOnly.isDBSNP());
        assertTrue(hasCOSMICOnly.isCOSMIC());

        final String none = "0 \t 1 \t 2 \t 3 \t 4 \t 5 \t 6 \t 7 \t 8 \t 9";
        final SomaticVariant hasNone = SomaticVariantFactory.fromVCFLine(none);
        assertFalse(hasNone.isDBSNP());
        assertFalse(hasNone.isCOSMIC());
    }

    @Test
    public void correctMissense() {
        final String missense = "0 \t 1 \t 2 \t 3 \t 4 \t 5 \t 6 \t missense_variant \t 8 \t 9";
        final SomaticVariant isMissense = SomaticVariantFactory.fromVCFLine(missense);
        assertTrue(isMissense.hasConsequence(VariantConsequence.MISSENSE_VARIANT));

        final String notMissense = "0 \t 1 \t 2 \t 3 \t 4 \t 5 \t 6 \t intron_variant \t 8 \t 9";
        final SomaticVariant isNotMissense = SomaticVariantFactory.fromVCFLine(notMissense);
        assertFalse(isNotMissense.hasConsequence(VariantConsequence.MISSENSE_VARIANT));
    }
}