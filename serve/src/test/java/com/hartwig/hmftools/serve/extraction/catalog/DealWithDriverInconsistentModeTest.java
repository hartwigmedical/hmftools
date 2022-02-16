package com.hartwig.hmftools.serve.extraction.catalog;

import org.junit.Test;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

public class DealWithDriverInconsistentModeTest {

    @Test
    public void canTestFilterOnInconsistenties() {
        assertFalse(DealWithDriverInconsistentMode.filterOnInconsistenties(DealWithDriverInconsistentModeAnnotation.IGNORE));
        assertTrue(DealWithDriverInconsistentMode.filterOnInconsistenties(DealWithDriverInconsistentModeAnnotation.WARN_ONLY));
        assertTrue(DealWithDriverInconsistentMode.filterOnInconsistenties(DealWithDriverInconsistentModeAnnotation.FILTER));
    }
}