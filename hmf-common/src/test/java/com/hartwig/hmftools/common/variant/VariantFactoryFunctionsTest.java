package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import org.junit.Test;

public class VariantFactoryFunctionsTest {

    @Test
    public void canSplitSampleDataField() {
        final String typicalGermlineSampleData = "0/1:32,61:93:99:2156,0,1092";
        final String typicalSomaticSampleData = "0/1:48,15:64";
        final String truthSetSomaticSampleData = "0/1:.:567:415";

        assertEquals(5, VariantFactoryFunctions.splitSampleDataFields(typicalGermlineSampleData).length);
        assertEquals(3, VariantFactoryFunctions.splitSampleDataFields(typicalSomaticSampleData).length);
        assertEquals(4, VariantFactoryFunctions.splitSampleDataFields(truthSetSomaticSampleData).length);
    }

    @Test
    public void canDetermineAlleleFrequencies() {
        final AlleleFrequencyData simple = VariantFactoryFunctions.determineAlleleFrequencies("60,30");
        assertNotNull(simple);
        assertEquals(30, simple.alleleReadCount());
        assertEquals(90, simple.totalReadCount());

        final AlleleFrequencyData complex = VariantFactoryFunctions.determineAlleleFrequencies("60,30,10");
        assertNotNull(complex);
        assertEquals(30, complex.alleleReadCount());
        assertEquals(100, complex.totalReadCount());

        final AlleleFrequencyData invalid = VariantFactoryFunctions.determineAlleleFrequencies("60");
        assertNull(invalid);
    }
}