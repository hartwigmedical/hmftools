package com.hartwig.hmftools.serve.actionability;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.hartwig.hmftools.serve.util.RefGenomeVersion;

import org.junit.Test;

public class ActionableEventsLoaderTest {

    @Test
    public void canLoadFromTestDir() throws IOException {
        assertNotNull(ActionableEventsLoader.readFromDir(ActionabilityTestUtil.TEST_ACTIONABILITY_DIR, RefGenomeVersion.V37));
    }
}