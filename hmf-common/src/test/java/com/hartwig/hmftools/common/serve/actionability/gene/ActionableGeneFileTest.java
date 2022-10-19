package com.hartwig.hmftools.common.serve.actionability.gene;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.serve.actionability.ActionabilityTestUtil;

import org.junit.Test;

public class ActionableGeneFileTest {

    @Test
    public void canReadFromFileAndConvert() throws IOException {
        String actionableGeneTsv =
                ActionableGeneFile.actionableGeneTsvPath(ActionabilityTestUtil.TEST_SERVE_OUTPUT_DIR, RefGenomeVersion.V37);
        List<ActionableGene> actionableGenes = ActionableGeneFile.read(actionableGeneTsv);

        assertEquals(7, actionableGenes.size());

        List<String> lines = ActionableGeneFile.toLines(actionableGenes);
        List<ActionableGene> regeneratedGenes = ActionableGeneFile.fromLines(lines);
        List<String> regeneratedLines = ActionableGeneFile.toLines(regeneratedGenes);
        assertEquals(lines.size(), regeneratedLines.size());

        for (int i = 0; i < lines.size(); i++) {
            assertEquals(lines.get(i), regeneratedLines.get(i));
        }
    }
}