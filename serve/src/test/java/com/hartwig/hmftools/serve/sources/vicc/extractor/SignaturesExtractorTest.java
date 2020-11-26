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

import org.junit.Ignore;
import org.junit.Test;

public class SignaturesExtractorTest {

    @Test
    @Ignore
    public void canExtractSignatureNameUnknown() {
        SignaturesExtractor signaturesExtractor = new SignaturesExtractor();
        ViccEntry viccEntry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.CIVIC,
                "any",
                "Other Biomarkers","description",
                "Microsatellite Instability-High", "chromosome", "pos");
        Feature feature = viccEntry.features().get(0);
        SignatureName signatureName = SignaturesExtractor.extractSignatureName(feature
                .name());

        //extract signature name
        assertEquals(SignatureName.MICROSATELLITE_UNSTABLE, signatureName);

        Map<Feature, SignatureName> signaturesPerFeature = Maps.newHashMap();
        signaturesPerFeature.put(feature, signatureName);

        assertEquals(signaturesPerFeature, signaturesExtractor.extractSignatures(viccEntry));
    }

    @Test
    public void canExtractSignatureUnknown() {
        SignaturesExtractor signaturesExtractor = new SignaturesExtractor();
        ViccEntry viccEntry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.CIVIC,
                "any",
                "Other Biomarkers","description",
                "Tum", "chromosome", "pos");
        Feature feature = viccEntry.features().get(0);
        SignatureName signatureName = SignaturesExtractor.extractSignatureName(feature
                .name());

        //extract signature name
        assertNull(signatureName);

        Map<Feature, SignatureName> signaturesPerFeature = Maps.newHashMap();

        assertEquals(signaturesPerFeature, signaturesExtractor.extractSignatures(viccEntry));
    }
}