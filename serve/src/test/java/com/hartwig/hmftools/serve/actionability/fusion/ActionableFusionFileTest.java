package com.hartwig.hmftools.serve.actionability.fusion;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.serve.actionability.ActionabilityTestUtil;

import org.junit.Test;

public class ActionableFusionFileTest {

    @Test
    public void canLoadActionableFusions() throws IOException {
        String actionableFusionTsv = ActionableFusionFile.actionableFusionTsvFilePath(ActionabilityTestUtil.SERVE_ACTIONABILITY_DIR);
        List<ActionableFusion> actionableFusions = ActionableFusionFile.loadFromActionableFusionTsv(actionableFusionTsv);

        assertEquals(2, actionableFusions.size());
    }
}