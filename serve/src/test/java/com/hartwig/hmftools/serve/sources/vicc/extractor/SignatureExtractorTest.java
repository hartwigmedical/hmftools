package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.actionability.signature.SignatureName;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class SignatureExtractorTest {

    private static final String MSI = "msi";
    private static final String TMB = "tmb";
    private static final String HRD = "hrd";

    @Test
    public void canExtractMicrosatelliteUnstableSignature() {
        SignatureExtractor signatureExtractor = buildTestExtractor();
        SignatureName signature = signatureExtractor.extract(EventType.SIGNATURE, MSI);

        assertNotNull(signature);
        assertEquals(SignatureName.MICROSATELLITE_UNSTABLE, signature);
    }

    @Test
    public void canExtractHighTumorMutationalLoadSignature() {
        SignatureExtractor signatureExtractor = buildTestExtractor();
        SignatureName signature = signatureExtractor.extract(EventType.SIGNATURE, TMB);

        assertNotNull(signature);
        assertEquals(SignatureName.HIGH_TUMOR_MUTATIONAL_LOAD, signature);
    }

    @Test
    public void canExtractHrDeficientSignature() {
        SignatureExtractor signatureExtractor = buildTestExtractor();
        SignatureName signature = signatureExtractor.extract(EventType.SIGNATURE, HRD);

        assertNotNull(signature);
        assertEquals(SignatureName.HOMOLOGOUS_RECOMBINATION_DEFICIENCY, signature);
    }

    @Test
    public void canFilterUnknownSignature() {
        SignatureExtractor signatureExtractor = buildTestExtractor();

        assertNull(signatureExtractor.extract(EventType.SIGNATURE, "Not a signature"));
    }

    @Test
    public void canFilterWrongTypes() {
        SignatureExtractor signatureExtractor = buildTestExtractor();

        assertNull(signatureExtractor.extract(EventType.COMPLEX, MSI));
    }

    @NotNull
    private static SignatureExtractor buildTestExtractor() {
        return new SignatureExtractor(Sets.newHashSet(MSI), Sets.newHashSet(TMB), Sets.newHashSet(HRD));
    }
}