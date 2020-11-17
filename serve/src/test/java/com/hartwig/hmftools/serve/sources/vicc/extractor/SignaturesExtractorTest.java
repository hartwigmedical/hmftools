package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.serve.actionability.signature.SignatureName;

import org.junit.Test;

public class SignaturesExtractorTest {

    @Test
    public void canExtractSignatureName() {
        assertEquals(SignatureName.MICROSATELLITE_UNSTABLE, SignaturesExtractor.extractSignatureName("Microsatellite Instability-High"));
        assertNull(SignaturesExtractor.extractSignatureName("abc"));
    }
}