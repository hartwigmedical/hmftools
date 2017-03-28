package com.hartwig.hmftools.healthchecker.context;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.FileNotFoundException;

import com.google.common.io.Resources;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class MetaDataResolverTest {

    private static final String RESOURCE_DIR = Resources.getResource("MetaDataResolver").getPath();

    @Test
    public void noMetaDataReturnsNull() throws FileNotFoundException {
        final String noMetaDataRunDir = RESOURCE_DIR + File.separator + "RunDirNoMetaData";
        assertNull(MetaDataResolver.fromMetaDataFile(noMetaDataRunDir));
    }

    @Test
    public void canResolveSingleSampleMetaData() throws FileNotFoundException {
        final String runName = "RunDirSingleSample";
        final String singleSampleRunDir = RESOURCE_DIR + File.separator + runName;
        final RunContext runContext = MetaDataResolver.fromMetaDataFile(singleSampleRunDir);

        assertNotNull(runContext);
        assertFalse(runContext.isSomaticRun());
        assertEquals("GIAB", runContext.refSample());
        assertEquals(Strings.EMPTY, runContext.tumorSample());
        assertEquals(runName, runContext.setName());
        assertEquals(singleSampleRunDir, runContext.runDirectory());
    }

    @Test
    public void canResolveSomaticMetaData() throws FileNotFoundException {
        final String runName = "RunDirSomatic";
        final String somaticRunDir = RESOURCE_DIR + File.separator + runName;
        final RunContext runContext = MetaDataResolver.fromMetaDataFile(somaticRunDir);

        assertNotNull(runContext);
        assertTrue(runContext.isSomaticRun());
        assertEquals("CPCT12345678R", runContext.refSample());
        assertEquals("CPCT12345678T", runContext.tumorSample());
        assertEquals(runName, runContext.setName());
        assertEquals(somaticRunDir, runContext.runDirectory());
    }
}
