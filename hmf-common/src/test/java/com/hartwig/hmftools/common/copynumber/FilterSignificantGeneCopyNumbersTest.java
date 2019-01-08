package com.hartwig.hmftools.common.copynumber;

import static org.junit.Assert.*;

import org.junit.Test;

public class FilterSignificantGeneCopyNumbersTest {

    @Test
    public void canDetermineSignificantEvent() {
        assertTrue(FilterSignificantGeneCopyNumbers.isSignificant(1, 0));
        assertTrue(FilterSignificantGeneCopyNumbers.isSignificant(2, 0));
        assertFalse(FilterSignificantGeneCopyNumbers.isSignificant(2, 1));

        assertFalse(FilterSignificantGeneCopyNumbers.isSignificant(2, 2));
        assertTrue(FilterSignificantGeneCopyNumbers.isSignificant(2, 20));
        assertFalse(FilterSignificantGeneCopyNumbers.isSignificant(10, 20));
    }

}