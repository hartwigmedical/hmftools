package com.hartwig.hmftools.serve.actionability.variant;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.serve.actionability.ActionabilityTestUtil;

import org.junit.Test;

public class ActionableVariantFileTest {

    @Test
    public void canReadFromFileAndConvert() throws IOException {
        String actionableVariantTsv = ActionableVariantFile.actionableVariantTsvPath(ActionabilityTestUtil.SERVE_ACTIONABILITY_DIR);
        List<ActionableVariant> actionableVariants = ActionableVariantFile.loadFromActionableVariantTsv(actionableVariantTsv);

        assertEquals(1, actionableVariants.size());

        List<String> lines = ActionableVariantFile.toLines(actionableVariants);
        List<ActionableVariant> regeneratedVariants = ActionableVariantFile.fromLines(lines);
        List<String> regeneratedLines = ActionableVariantFile.toLines(regeneratedVariants);
        assertEquals(lines.size(), regeneratedLines.size());

        for (int i = 0; i < lines.size(); i++) {
            assertEquals(lines.get(i), regeneratedLines.get(i));
        }
    }
}