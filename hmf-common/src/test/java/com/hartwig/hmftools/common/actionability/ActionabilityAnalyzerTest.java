package com.hartwig.hmftools.common.actionability;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class ActionabilityAnalyzerTest {

    private static final String KNOWLEDGEBASE_PATH = Resources.getResource("actionability").getPath();

    @Test
    public void canGenerateFromTestKnowledgebase() throws IOException {
        ActionabilityAnalyzer actionabilityAnalyzer = ActionabilityAnalyzer.fromKnowledgebase(KNOWLEDGEBASE_PATH);

        assertEquals(2, actionabilityAnalyzer.variantAnalyzer().actionableGenes().size());
        assertEquals(1, actionabilityAnalyzer.cnvAnalyzer().actionableGenes().size());
        assertEquals(4, actionabilityAnalyzer.fusionAnalyzer().actionableGenes().size());
    }
}