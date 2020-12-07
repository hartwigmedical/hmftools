package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.Map;

import com.hartwig.hmftools.serve.actionability.signature.SignatureName;
import com.hartwig.hmftools.serve.sources.vicc.ViccTestFactory;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.junit.Test;

public class SignatureExtractorTest {

    @Test
    public void canExtractMSISignature() {
        SignatureExtractor signatureExtractor = new SignatureExtractor();
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("Other Biomarkers", "Microsatellite Instability-High");

        Map<Feature, SignatureName> signaturesPerFeature = signatureExtractor.extract(viccEntry);

        assertEquals(1, signaturesPerFeature.size());
        assertEquals(SignatureName.MICROSATELLITE_UNSTABLE, signaturesPerFeature.get(viccEntry.features().get(0)));
    }

    @Test
    public void doesNotFailOnUnknownSignature() {
        SignatureExtractor signatureExtractor = new SignatureExtractor();
        ViccEntry viccEntry = ViccTestFactory.testEntryWithGeneAndEvent("Other Biomarkers", "Not a signature");

        Map<Feature, SignatureName> signaturesPerFeature = signatureExtractor.extract(viccEntry);

        assertTrue(signaturesPerFeature.isEmpty());
    }

    @Test
    public void canExtractSignatureName() {
        assertEquals(SignatureName.MICROSATELLITE_UNSTABLE, SignatureExtractor.extractSignatureName("Microsatellite Instability-High"));

        assertNull(SignatureExtractor.extractSignatureName("abc"));
    }
}