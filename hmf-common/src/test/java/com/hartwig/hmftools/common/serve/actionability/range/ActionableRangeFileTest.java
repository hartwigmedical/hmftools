package com.hartwig.hmftools.common.serve.actionability.range;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.serve.actionability.ActionabilityTestUtil;

import org.junit.Test;

public class ActionableRangeFileTest {

    @Test
    public void canReadFromFileAndConvert() throws IOException {
        String actionableRangeTsv =
                ActionableRangeFile.actionableRangeTsvPath(ActionabilityTestUtil.TEST_SERVE_OUTPUT_DIR, RefGenomeVersion.V37);
        List<ActionableRange> actionableRanges = ActionableRangeFile.read(actionableRangeTsv);

        assertEquals(2, actionableRanges.size());

        List<String> lines = ActionableRangeFile.toLines(actionableRanges);
        List<ActionableRange> regeneratedRanges = ActionableRangeFile.fromLines(lines);
        List<String> regeneratedLines = ActionableRangeFile.toLines(regeneratedRanges);
        assertEquals(lines.size(), regeneratedLines.size());

        for (int i = 0; i < lines.size(); i++) {
            assertEquals(lines.get(i), regeneratedLines.get(i));
        }
    }
}