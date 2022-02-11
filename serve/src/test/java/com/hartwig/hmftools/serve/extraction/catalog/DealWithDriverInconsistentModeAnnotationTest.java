package com.hartwig.hmftools.serve.extraction.catalog;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class DealWithDriverInconsistentModeAnnotationTest {

    @Test
    public void canExtractAnnotation() {
        DealWithDriverInconsistentModeAnnotation annotationFilter =
                DealWithDriverInconsistentModeAnnotation.extractDealWithDriverInconsistentMode("filter");
        assertEquals(DealWithDriverInconsistentModeAnnotation.FILTER, annotationFilter);
        assertTrue(DealWithDriverInconsistentModeAnnotation.FILTER.logging());

        DealWithDriverInconsistentModeAnnotation annotationIgnore =
                DealWithDriverInconsistentModeAnnotation.extractDealWithDriverInconsistentMode("ignore");
        assertEquals(DealWithDriverInconsistentModeAnnotation.IGNORE, annotationIgnore);
        assertFalse(DealWithDriverInconsistentModeAnnotation.IGNORE.logging());

        DealWithDriverInconsistentModeAnnotation annotationWarnOnly =
                DealWithDriverInconsistentModeAnnotation.extractDealWithDriverInconsistentMode("warn_only");
        assertEquals(DealWithDriverInconsistentModeAnnotation.WARN_ONLY, annotationWarnOnly);
        assertTrue(DealWithDriverInconsistentModeAnnotation.WARN_ONLY.logging());

    }

    @Test(expected = IllegalStateException.class)
    public void hasUnknownInput() {
        DealWithDriverInconsistentModeAnnotation.extractDealWithDriverInconsistentMode("AB");
    }
}