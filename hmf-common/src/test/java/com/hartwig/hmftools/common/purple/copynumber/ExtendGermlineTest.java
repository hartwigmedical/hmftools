package com.hartwig.hmftools.common.purple.copynumber;

import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertFalse;

import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ExtendGermlineTest {

    @Test
    public void testValidAmplification() {
        assertTrue(ExtendGermline.isValidAmplification(4, 5, create(SegmentSupport.BND), null));
        assertTrue(ExtendGermline.isValidAmplification(4, 5, create(SegmentSupport.INS), create(SegmentSupport.NONE)));
        assertTrue(ExtendGermline.isValidAmplification(4, 5, create(SegmentSupport.INV), create(SegmentSupport.DEL)));
        assertTrue(ExtendGermline.isValidAmplification(4, 5, create(SegmentSupport.NONE), create(SegmentSupport.DUP)));
        assertTrue(ExtendGermline.isValidAmplification(4.1, 5.1, create(SegmentSupport.NONE), create(SegmentSupport.MULTIPLE)));

        assertFalse(ExtendGermline.isValidAmplification(1, 5, create(SegmentSupport.BND), create(SegmentSupport.CENTROMERE)));
        assertFalse(ExtendGermline.isValidAmplification(1, 5, create(SegmentSupport.CENTROMERE), create(SegmentSupport.BND)));
        assertFalse(ExtendGermline.isValidAmplification(1, 5, create(SegmentSupport.NONE), create(SegmentSupport.NONE)));

        assertFalse(ExtendGermline.isValidAmplification(1, 4.9, create(SegmentSupport.BND), create(SegmentSupport.NONE)));
        assertFalse(ExtendGermline.isValidAmplification(4.1, 5, create(SegmentSupport.BND), create(SegmentSupport.NONE)));
    }

    @NotNull
    private static FittedRegion create(@NotNull final SegmentSupport support) {
        return PurpleDatamodelTest.createDefaultFittedRegion("1", 1, 1000).support(support).build();
    }
}
