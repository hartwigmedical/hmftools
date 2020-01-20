package com.hartwig.hmftools.protect.actionability.cnv;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class SignificantGeneCopyNumberFilterTest {

    @Test
    public void canDetermineSignificantEvent() {
        assertTrue(SignificantGeneCopyNumberFilter.isSignificant(1, 0));
        assertTrue(SignificantGeneCopyNumberFilter.isSignificant(2, 0.3));
        assertFalse(SignificantGeneCopyNumberFilter.isSignificant(2, 1));

        assertFalse(SignificantGeneCopyNumberFilter.isSignificant(2, 2));
        assertTrue(SignificantGeneCopyNumberFilter.isSignificant(2, 20));
        assertFalse(SignificantGeneCopyNumberFilter.isSignificant(10, 20));
    }
}