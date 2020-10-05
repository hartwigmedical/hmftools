package com.hartwig.hmftools.serve.actionability.signature;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.serve.actionability.ActionabilityTestUtil;

import org.junit.Test;

public class ActionableSignatureFileTest {

    @Test
    public void canReadFromFileAndConvert() throws IOException {
        String actionableSignatureTsv = ActionableSignatureFile.actionableSignatureTsvPath(ActionabilityTestUtil.SERVE_ACTIONABILITY_DIR);
        List<ActionableSignature> actionableSignatures = ActionableSignatureFile.read(actionableSignatureTsv);

        assertEquals(1, actionableSignatures.size());

        List<String> lines = ActionableSignatureFile.toLines(actionableSignatures);
        List<ActionableSignature> regeneratedSignatures = ActionableSignatureFile.fromLines(lines);
        List<String> regeneratedLines = ActionableSignatureFile.toLines(regeneratedSignatures);
        assertEquals(lines.size(), regeneratedLines.size());

        for (int i = 0; i < lines.size(); i++) {
            assertEquals(lines.get(i), regeneratedLines.get(i));
        }
    }
}