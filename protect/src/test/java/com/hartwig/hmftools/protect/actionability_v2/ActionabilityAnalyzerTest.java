package com.hartwig.hmftools.protect.actionability_v2;

import static org.junit.Assert.*;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class ActionabilityAnalyzerTest {

    private static final String KNOWLEDGEBASE_DIRECTORY_V2 = Resources.getResource("actionabilty_v2").getPath();

    @Test
    public void canGenerateFromTestKnowledgebase() throws IOException {
         ActionabilityAnalyzer actionabilityAnalyzer = ActionabilityAnalyzer.fromKnowledgebase(KNOWLEDGEBASE_DIRECTORY_V2);

    }

}