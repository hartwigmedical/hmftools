package com.hartwig.hmftools.serve.actionability;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class ActionabilityAnalyzerTest {

    private static final String KNOWLEDGEBASE_DIRECTORY_V2 = Resources.getResource("actionability").getPath();

    @Test
    public void canGenerateFromTestKnowledgebase() throws IOException {
         ActionabilityAnalyzer.fromKnowledgebase(KNOWLEDGEBASE_DIRECTORY_V2);
    }
}