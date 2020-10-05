package com.hartwig.hmftools.serve.actionability.gene;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.serve.actionability.ActionabilityTestUtil;

import org.junit.Test;

public class ActionableGeneFileTest {

    @Test
    public void canLoadActionableGenes() throws IOException {
        String actionableGeneTsv = ActionableGeneFile.actionableGeneTsvFilePath(ActionabilityTestUtil.SERVE_ACTIONABILITY_DIR);
        List<ActionableGene> actionableGenes = ActionableGeneFile.loadFromActionableGeneTsv(actionableGeneTsv);

        assertEquals(7, actionableGenes.size());
    }
}