package com.hartwig.hmftools.serve.actionability.variant;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.serve.actionability.ActionabilityTestUtil;

import org.junit.Test;

public class VariantEvidenceAnalyzerFactoryTest {

    @Test
    public void canLoadActionableVariants() throws IOException {
        String actionableVariantTsv =
                VariantEvidenceAnalyzerFactory.actionableVariantTsvFilePath(ActionabilityTestUtil.SERVE_ACTIONABILITY_DIR);
        List<ActionableVariant> actionableVariants = VariantEvidenceAnalyzerFactory.loadFromActionableVariantTsv(actionableVariantTsv);

        assertEquals(1, actionableVariants.size());
    }
}