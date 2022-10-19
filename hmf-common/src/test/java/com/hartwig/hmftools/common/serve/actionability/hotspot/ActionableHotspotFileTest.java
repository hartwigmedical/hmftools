package com.hartwig.hmftools.common.serve.actionability.hotspot;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.serve.actionability.ActionabilityTestUtil;

import org.junit.Test;

public class ActionableHotspotFileTest {

    @Test
    public void canReadFromFileAndConvert() throws IOException {
        String actionableHotspotTsv =
                ActionableHotspotFile.actionableHotspotTsvPath(ActionabilityTestUtil.TEST_SERVE_OUTPUT_DIR, RefGenomeVersion.V37);
        List<ActionableHotspot> actionableHotspots = ActionableHotspotFile.read(actionableHotspotTsv);

        assertEquals(2, actionableHotspots.size());

        List<String> lines = ActionableHotspotFile.toLines(actionableHotspots);
        List<ActionableHotspot> regeneratedHotspots = ActionableHotspotFile.fromLines(lines);
        List<String> regeneratedLines = ActionableHotspotFile.toLines(regeneratedHotspots);
        assertEquals(lines.size(), regeneratedLines.size());

        for (int i = 0; i < lines.size(); i++) {
            assertEquals(lines.get(i), regeneratedLines.get(i));
        }
    }
}