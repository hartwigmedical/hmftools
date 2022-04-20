package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class CanonicalAnnotationTest {

    @Test
    public void testTrimEnsembleTranscriptId() {
        assertEquals("ENST00000361570", CanonicalAnnotation.trimEnsembleVersion("ENST00000361570"));
        assertEquals("ENST00000361570", CanonicalAnnotation.trimEnsembleVersion("ENST00000361570.v8"));
    }
}
