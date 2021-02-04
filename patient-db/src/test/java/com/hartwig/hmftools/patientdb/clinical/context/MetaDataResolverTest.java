package com.hartwig.hmftools.patientdb.clinical.context;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class MetaDataResolverTest {

    private static final String RESOURCE_DIR = Resources.getResource("context").getPath();

    @Test
    public void noMetaDataReturnsNull() throws IOException {
        String noMetaDataRunDir = RESOURCE_DIR + File.separator + "RunDirNoMetaData";
        assertNull(MetaDataResolver.fromMetaDataFile(noMetaDataRunDir, ""));
    }

    @Test
    public void noRefSampleReturnsNull() throws IOException{
        String noRefSampleRunDir = RESOURCE_DIR + File.separator + "RunDirNoRefSample";
        assertNull(MetaDataResolver.fromMetaDataFile(noRefSampleRunDir, ""));
    }

    @Test
    public void noSetNameReturnsNull() throws IOException{
        String noSetNameRunDir = RESOURCE_DIR + File.separator + "RunDirNoSetName";
        assertNull(MetaDataResolver.fromMetaDataFile(noSetNameRunDir, ""));
    }

    @Test
    public void canResolveMetaDataFilePV5() throws IOException{
        String noSetNameRunDir = RESOURCE_DIR + File.separator + "RunDirP5";
        assertNotNull(MetaDataResolver.fromMetaDataFile(noSetNameRunDir, "pipeline.version"));
    }

    @Test
    public void canResolveSomaticMetaDataP4() throws IOException{
        String setName = "RunDirSomatic";
        String runDirectory = RESOURCE_DIR + File.separator + setName;
        RunContext runContext = MetaDataResolver.fromMetaDataFile(runDirectory, "");

        assertNotNull(runContext);
        assertEquals("CPCT12345678R", runContext.refSample());
        assertEquals("CPCT12345678T", runContext.tumorSample());
        assertEquals(setName, runContext.setName());
        assertEquals(runDirectory, runContext.runDirectory());
        assertEquals(Strings.EMPTY, runContext.tumorBarcodeSample());
    }

    @Test
    public void canResolveSomaticMetaDataP5() throws IOException{
        String setName = "RunDirP5";
        String runDirectory = RESOURCE_DIR + File.separator + setName;
        RunContext runContext = MetaDataResolver.fromMetaDataFile(runDirectory, "pipeline.version");

        assertNotNull(runContext);
        assertEquals("CPCT12345678R", runContext.refSample());
        assertEquals("CPCT12345678T", runContext.tumorSample());
        assertEquals(setName, runContext.setName());
        assertEquals(runDirectory, runContext.runDirectory());
        assertEquals("AB456", runContext.tumorBarcodeSample());
    }

    @Test
    public void canResolveSomaticMetaDataPostP5_15() throws IOException{
        String setName = "RunDirPostP5_15";
        String runDirectory = RESOURCE_DIR + File.separator + setName;
        RunContext runContext = MetaDataResolver.fromMetaDataFile(runDirectory, "pipeline.version");

        assertNotNull(runContext);
        assertEquals("CPCT12345678R", runContext.refSample());
        assertEquals("CPCT12345678T", runContext.tumorSample());
        assertEquals(setName, runContext.setName());
        assertEquals(runDirectory, runContext.runDirectory());
        assertEquals("AB456", runContext.tumorBarcodeSample());
    }
}
