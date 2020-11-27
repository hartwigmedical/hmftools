package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.serve.actionability.signature.SignatureName;
import com.hartwig.hmftools.serve.sources.vicc.ViccTestFactory;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.junit.Test;

public class SignaturesExtractorTest {

    @Test
    public void canExtractSignatureName() {
        assertEquals(SignatureName.MICROSATELLITE_UNSTABLE, SignaturesExtractor.extractSignatureName("Microsatellite Instability-High"));
        assertNull(SignaturesExtractor.extractSignatureName("abc"));
    }

    @Test
    public void canExtractSignatureNameUnknown() {
        SignaturesExtractor signaturesExtractor = new SignaturesExtractor();
        ViccEntry viccEntry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.CIVIC,
                "any",
                "Other Biomarkers",
                "Microsatellite Instability-High",
                "chromosome",
                "pos",
                null);
        Feature feature = viccEntry.features().get(0);
        SignatureName signatureName = SignaturesExtractor.extractSignatureName(feature.name());

        Map<Feature, SignatureName> signaturesPerFeature = Maps.newHashMap();
        signaturesPerFeature.put(feature, signatureName);

        assertEquals(signaturesPerFeature, signaturesExtractor.extractSignatures(viccEntry));
    }

    @Test
    public void canExtractSignatureUnknown() {
        SignaturesExtractor signaturesExtractor = new SignaturesExtractor();
        ViccEntry viccEntry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.CIVIC,
                "any",
                "Other Biomarkers",
                "Tum",
                "chromosome",
                "pos",
                null);

        Map<Feature, SignatureName> signaturesPerFeature = Maps.newHashMap();

        assertEquals(signaturesPerFeature, signaturesExtractor.extractSignatures(viccEntry));
    }
}