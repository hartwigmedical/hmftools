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

        assertEquals(5, VariantFactoryFunctions.splitSampleDataFields(typicalGermlineSampleData).length);
        assertEquals(3, VariantFactoryFunctions.splitSampleDataFields(typicalSomaticSampleData).length);
    }

    @Test
    public void canDetermineAlleleFrequencies() {
        final ReadCount simpleReadCount = VariantFactoryFunctions.analyzeAlleleFrequencies("60,30");
        assertNotNull(simpleReadCount);
        assertEquals(30, simpleReadCount.alleleReadCount());
        assertEquals(90, simpleReadCount.totalReadCount());

        final ReadCount complexReadCount = VariantFactoryFunctions.analyzeAlleleFrequencies("60,30,10");
        assertNotNull(complexReadCount);
        assertEquals(30, complexReadCount.alleleReadCount());
        assertEquals(100, complexReadCount.totalReadCount());

        final ReadCount invalidReadCount = VariantFactoryFunctions.analyzeAlleleFrequencies("60");
        assertNull(invalidReadCount);
    }
}