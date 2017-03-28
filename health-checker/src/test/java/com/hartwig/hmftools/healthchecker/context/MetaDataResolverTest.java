package com.hartwig.hmftools.healthchecker.context;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.File;

import com.google.common.io.Resources;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class MetaDataResolverTest {

    private static final String RESOURCE_DIR = Resources.getResource("MetaDataResolver").getPath();

    @Test
    public void noMetaDataReturnsNull() {
        final String noMetaDataRunDir = RESOURCE_DIR + File.separator + "RunDirNoMetaData";
        assertNull(MetaDataResolver.fromMetaDataFile(noMetaDataRunDir));
    }

    @Test
    public void noRefSampleReturnsNull() {
        final String noRefSampleRunDir = RESOURCE_DIR + File.separator + "RunDirNoRefSample";
        assertNull(MetaDataResolver.fromMetaDataFile(noRefSampleRunDir));
    }

    @Test
    public void noSetNameReturnsNull() {
        final String noSetNameRunDir = RESOURCE_DIR + File.separator + "RunDirNoSetName";
        assertNull(MetaDataResolver.fromMetaDataFile(noSetNameRunDir));
    }

    @Test
    public void canResolveSingleSampleMetaDataWithExplicitNoTumorSample() {
        testSingleSample("RunDirSingleSampleWithTumorSample");
    }

    @Test
    public void canResolveSingleSampleMetaDataWithoutTumorSample() {
        testSingleSample("RunDirSingleSampleNoTumorSample");
    }

    private static void testSingleSample(@NotNull final String setName) {
        final String runDirectory = RESOURCE_DIR + File.separator + setName;
        final RunContext runContext = MetaDataResolver.fromMetaDataFile(runDirectory);

        assertNotNull(runContext);
        assertFalse(runContext.isSomaticRun());
        assertEquals("GIAB", runContext.refSample());
        assertEquals(Strings.EMPTY, runContext.tumorSample());
        assertEquals(setName, runContext.setName());
        assertEquals(runDirectory, runContext.runDirectory());
    }

    @Test
    public void canResolveSomaticMetaData() {
        final String setName = "RunDirSomatic";
        final String runDirectory = RESOURCE_DIR + File.separator + setName;
        final RunContext runContext = MetaDataResolver.fromMetaDataFile(runDirectory);

        assertNotNull(runContext);
        assertTrue(runContext.isSomaticRun());
        assertEquals("CPCT12345678R", runContext.refSample());
        assertEquals("CPCT12345678T", runContext.tumorSample());
        assertEquals(setName, runContext.setName());
        assertEquals(runDirectory, runContext.runDirectory());
    }
}
