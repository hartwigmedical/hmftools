package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class GermlineVariantFactoryTest {

    private static final String VALID_SAMPLE_DATA = "0/1:20,10:35:whatever";
    private static final String INVALID_SAMPLE_DATA = "0/1:20,10";

    @Test
    public void extractFromNormalVCFLine() {
        final String somewhatValidVCFLine =
                "1 \t 2 \t 3 \t 4 \t 5 \t 6 \t 7 \t 8 \t 9 \t " + VALID_SAMPLE_DATA + " \t " + VALID_SAMPLE_DATA;
        GermlineVariant variant = GermlineVariantFactory.fromVCFLine(somewhatValidVCFLine);
        assertNotNull(variant);
        assertNotNull(variant.refData());
        assertNotNull(variant.tumorData());
    }

    @Test
    public void extractFromInvalidVCFLine() {
        final String invalidVCFLine = "1 \t 2 \t 3 \t 4 \t 5 \t 6 \t 7 \t 8 \t 9 ";
        assertNull(GermlineVariantFactory.fromVCFLine(invalidVCFLine));
    }

    @Test
    public void extractFromInvalidSampleVCFLine() {
        final String invalidVCFLine =
                "1 \t 2 \t 3 \t 4 \t 5 \t 6 \t 7 \t 8 \t 9 \t " + INVALID_SAMPLE_DATA + " \t " + INVALID_SAMPLE_DATA;
        assertNull(GermlineVariantFactory.fromVCFLine(invalidVCFLine));
    }

    @Test
    public void canGenerateGermlineSampleData() {
        final GermlineSampleData validSampleData = GermlineVariantFactory.fromSampleData(VALID_SAMPLE_DATA);
        assertNotNull(validSampleData);
        assertEquals("0/1", validSampleData.genoType());
        assertEquals(10, validSampleData.alleleReadCount());
        assertEquals(30, validSampleData.totalReadCount());

        assertNull(GermlineVariantFactory.fromSampleData(Strings.EMPTY));
        assertNull(GermlineVariantFactory.fromSampleData("hi:there"));
    }
}