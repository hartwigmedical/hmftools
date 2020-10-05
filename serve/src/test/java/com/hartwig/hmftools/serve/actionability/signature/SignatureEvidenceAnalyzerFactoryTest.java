package com.hartwig.hmftools.serve.actionability.signature;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.serve.actionability.ActionabilityTestUtil;

import org.junit.Test;

public class SignatureEvidenceAnalyzerFactoryTest {

    @Test
    public void canLoadActionableSignatures() throws IOException {
        String actionableSignatureTsv =
                SignatureEvidenceAnalyzerFactory.actionableSignatureTsvFilePath(ActionabilityTestUtil.SERVE_ACTIONABILITY_DIR);
        List<ActionableSignature> actionableSignatures =
                SignatureEvidenceAnalyzerFactory.loadFromActionableSignatureTsv(actionableSignatureTsv);

        assertEquals(1, actionableSignatures.size());
    }
}