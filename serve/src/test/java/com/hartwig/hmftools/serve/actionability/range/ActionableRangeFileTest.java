package com.hartwig.hmftools.serve.actionability.range;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.serve.actionability.ActionabilityTestUtil;

import org.junit.Test;

public class ActionableRangeFileTest {

    @Test
    public void canLoadActionableRanges() throws IOException {
        String actionableRangeTsv = ActionableRangeFile.actionableRangeTsvFilePath(ActionabilityTestUtil.SERVE_ACTIONABILITY_DIR);
        List<ActionableRange> actionableRanges = ActionableRangeFile.loadFromActionableRangeTsv(actionableRangeTsv);

        assertEquals(2, actionableRanges.size());
    }
}