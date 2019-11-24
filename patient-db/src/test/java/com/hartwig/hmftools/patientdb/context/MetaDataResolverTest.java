package com.hartwig.hmftools.patientdb.context;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.io.File;

import com.google.common.io.Resources;

import org.junit.Test;

public class MetaDataResolverTest {

    private static final String RESOURCE_DIR = Resources.getResource("context").getPath();

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
    public void canResolveMetaDataFilePV5() {
        final String noSetNameRunDir = RESOURCE_DIR + File.separator + "RunDirP5";
        assertNotNull(MetaDataResolver.fromMetaDataFile(noSetNameRunDir));
    }

    @Test
    public void canResolveSomaticMetaDataP4() {
        final String setName = "RunDirSomatic";
        final String runDirectory = RESOURCE_DIR + File.separator + setName;
        final RunContext runContext = MetaDataResolver.fromMetaDataFile(runDirectory);

        assertNotNull(runContext);
        assertEquals("CPCT12345678R", runContext.refSample());
        assertEquals("CPCT12345678T", runContext.tumorSample());
        assertEquals(setName, runContext.setName());
        assertEquals(runDirectory, runContext.runDirectory());
        assertEquals("", runContext.tumorBarcodeSample());
    }

    @Test
    public void canResolveSomaticMetaDataP5() {
        final String setName = "RunDirP5";
        final String runDirectory = RESOURCE_DIR + File.separator + setName;
        final RunContext runContext = MetaDataResolver.fromMetaDataFile(runDirectory);

        assertNotNull(runContext);
        assertEquals("CPCT12345678R", runContext.refSample());
        assertEquals("CPCT12345678T", runContext.tumorSample());
        assertEquals(setName, runContext.setName());
        assertEquals(runDirectory, runContext.runDirectory());
        assertEquals("AB456", runContext.tumorBarcodeSample());
    }
}
