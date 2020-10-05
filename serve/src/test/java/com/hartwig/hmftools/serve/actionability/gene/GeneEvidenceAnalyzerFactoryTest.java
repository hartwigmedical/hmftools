package com.hartwig.hmftools.serve.actionability.gene;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.serve.actionability.ActionabilityTestUtil;

import org.junit.Test;

public class GeneEvidenceAnalyzerFactoryTest {

    @Test
    public void canLoadActionableGenes() throws IOException {
        String actionableGeneTsv =
                GeneEvidenceAnalyzerFactory.actionableGeneTsvFilePath(ActionabilityTestUtil.SERVE_ACTIONABILITY_DIR);
        List<ActionableGene> actionableGenes = GeneEvidenceAnalyzerFactory.loadFromActionableGeneTsv(actionableGeneTsv);

        assertEquals(7, actionableGenes.size());
    }
}