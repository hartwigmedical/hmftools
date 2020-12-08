package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.actionability.signature.SignatureName;

import org.junit.Test;

public class SignatureExtractorTest {

    @Test
    public void canExtractMSISignature() {
        SignatureExtractor signatureExtractor = new SignatureExtractor();
        SignatureName signature = signatureExtractor.extract(EventType.SIGNATURE, "Microsatellite Instability-High");

        assertNotNull(signature);
        assertEquals(SignatureName.MICROSATELLITE_UNSTABLE, signature);
    }

    @Test
    public void doesNotFailOnUnknownSignature() {
        SignatureExtractor signatureExtractor = new SignatureExtractor();
        SignatureName signature = signatureExtractor.extract(EventType.SIGNATURE, "Not a signature");

        assertNull(signature);
    }

    @Test
    public void canExtractSignatureName() {
        assertEquals(SignatureName.MICROSATELLITE_UNSTABLE, SignatureExtractor.extractSignatureName("Microsatellite Instability-High"));

        assertNull(SignatureExtractor.extractSignatureName("abc"));
    }
}