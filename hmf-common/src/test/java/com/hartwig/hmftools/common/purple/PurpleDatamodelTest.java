package com.hartwig.hmftools.common.purple;

import static org.junit.Assert.assertNotNull;

import org.junit.Test;

public class PurpleDatamodelTest {

    private static final String CHROMOSOME = "1";

    @Test
    public void testDefaultFittedRegion() {
        assertNotNull(PurpleTestUtils.createDefaultFittedRegion(CHROMOSOME, 1, 100).build());
    }

    @Test
    public void testDefaultCopyNumber() {
        assertNotNull(PurpleTestUtils.createCopyNumber(CHROMOSOME, 1, 100, 2).build());
    }

}
