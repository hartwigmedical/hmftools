package com.hartwig.hmftools.common.serve.actionability;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.junit.Test;

public class ActionableEventsLoaderTest {

    @Test
    public void canLoadFromTestDir() throws IOException {
        assertNotNull(ActionableEventsLoader.readFromDir(ActionabilityTestUtil.TEST_SERVE_OUTPUT_DIR, RefGenomeVersion.V37));
    }
}