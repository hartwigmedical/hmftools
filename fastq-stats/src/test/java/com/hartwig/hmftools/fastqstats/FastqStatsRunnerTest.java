package com.hartwig.hmftools.fastqstats;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class FastqStatsRunnerTest {

    private static final String RUN_DIR_PATH = Resources.getResource("170101_TEST_000_TESTFLOWCELL").getPath();

    @Test
    public void extractsCorrectFlowcellName() {
        final String flowcellName = FastqStatsRunner.getFlowcellName(RUN_DIR_PATH);
        assertEquals("TESTFLOWCELL", flowcellName);
    }

    @Test
    public void returnsUnknownFlowcellName() {
        final String flowcellName = FastqStatsRunner.getFlowcellName("whoops");
        assertEquals("Unknown", flowcellName);
    }

    @Test
    public void getsCorrectBaseCallsDir() throws IOException {
        final File baseCallsDir = FastqStatsRunner.getBaseCallsDir(RUN_DIR_PATH);
        assertTrue(baseCallsDir.exists());
        assertTrue(baseCallsDir.isDirectory());
        assertEquals("BaseCalls", baseCallsDir.getName());
        assertEquals("Intensities", baseCallsDir.getParentFile().getName());
        assertEquals("Data", baseCallsDir.getParentFile().getParentFile().getName());
    }

    @Test(expected = IOException.class)
    public void getBaseCallsDirThrowsForNonExistentDir() throws IOException {
        FastqStatsRunner.getBaseCallsDir("Whoops");
    }
}
