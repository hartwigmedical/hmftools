package com.hartwig.hmftools.knowledgebasegenerator.actionability;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class ActionabilityAnalyzerTest {

    private static final String KNOWLEDGEBASE_DIRECTORY_V2 = Resources.getResource("actionabilty").getPath();

    @Test
    public void canGenerateFromTestKnowledgebase() throws IOException {
         ActionabilityAnalyzer actionabilityAnalyzer = ActionabilityAnalyzer.fromKnowledgebase(KNOWLEDGEBASE_DIRECTORY_V2);

    }

}