package com.hartwig.hmftools.common.purple.copynumber;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;

import org.junit.Test;

public class PurpleCopyNumberTest {

    @Test
    public void testNegativeCopyNumber() {
        PurpleCopyNumber copyNumber = PurpleDatamodelTest.createCopyNumber("1", 1, 100, -12).build();
        assertEquals(0, copyNumber.value());
    }
}
