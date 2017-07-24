package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.GermlineVariantFactory.fromSampleData;
import static com.hartwig.hmftools.common.variant.GermlineVariantFactory.fromVCFLine;
import static com.hartwig.hmftools.common.variant.GermlineVariantFactory.sortReferenceAndTumor;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GermlineVariantFactoryTest {

    private static final String VALID_SAMPLE_DATA = "0/1:20,10:35:whatever";
    private static final String INVALID_SAMPLE_DATA = "0/1:20,10";

    @Test
    public void extractReferenceAndTumorSampleNames() {
        assertReferenceAndTumor("COLO829BL", "COLO829");
        assertReferenceAndTumor("SampleR", "SampleT");
        assertReferenceAndTumor("SampleR", "SampleTI");
        assertReferenceAndTumor("SampleR", "SampleTII");
    }

    private static void assertReferenceAndTumor(@NotNull final String expectedReference, @NotNull final String expectedTumor) {
        assertEquals(expectedReference, sortReferenceAndTumor(expectedReference, expectedTumor)[0]);
        assertEquals(expectedReference, sortReferenceAndTumor(expectedTumor, expectedReference)[0]);
        assertEquals(expectedTumor, sortReferenceAndTumor(expectedReference, expectedTumor)[1]);
        assertEquals(expectedTumor, sortReferenceAndTumor(expectedTumor, expectedReference)[1]);
    }

    @Test
    public void extractFromNormalVCFLine() {
        final String somewhatValidVCFLine =
                "1 \t 2 \t 3 \t 4 \t 5 \t 6 \t 7 \t 8 \t 9 \t " + VALID_SAMPLE_DATA + " \t " + VALID_SAMPLE_DATA;
        GermlineVariant variant = fromVCFLine(somewhatValidVCFLine);
        assertNotNull(variant);
        assertNotNull(variant.refData());
        assertNotNull(variant.tumorData());
    }

    @Test
    public void extractFromInvalidVCFLine() {
        final String invalidVCFLine = "1 \t 2 \t 3 \t 4 \t 5 \t 6 \t 7 \t 8 \t 9 ";
        assertNull(fromVCFLine(invalidVCFLine));
    }

    @Test
    public void extractFromInvalidSampleVCFLine() {
        final String invalidVCFLine = "1 \t 2 \t 3 \t 4 \t 5 \t 6 \t 7 \t 8 \t 9 \t " + INVALID_SAMPLE_DATA + " \t " + INVALID_SAMPLE_DATA;
        assertNull(fromVCFLine(invalidVCFLine));
    }

    @Test
    public void canGenerateGermlineSampleData() {
        final GermlineSampleData validSampleData = fromSampleData(VALID_SAMPLE_DATA);
        assertNotNull(validSampleData);
        assertEquals("0/1", validSampleData.genoType());
        assertEquals(10, validSampleData.alleleReadCount());
        assertEquals(30, validSampleData.totalReadCount());

        assertNull(fromSampleData(Strings.EMPTY));
        assertNull(fromSampleData("hi:there"));
    }
}