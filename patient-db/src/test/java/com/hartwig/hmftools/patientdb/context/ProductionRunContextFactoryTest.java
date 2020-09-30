package com.hartwig.hmftools.patientdb.context;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.utils.io.exception.MalformedFileException;

import org.junit.Test;

public class ProductionRunContextFactoryTest {

    private static final String RESOURCE_DIR = Resources.getResource("context").getPath();

    @Test
    public void picksMetadataWhenAvailable() throws IOException {
        String runDirectory = RESOURCE_DIR + File.separator + "RunDirSomatic";
        assertNotNull(ProductionRunContextFactory.fromRunDirectory(runDirectory));
    }

    @Test(expected = MalformedFileException.class)
    public void throwExceptionWhenNoMetaData() throws IOException {
        String runDirectory = RESOURCE_DIR + File.separator + "DoesNotExist";
        assertNotNull(ProductionRunContextFactory.fromRunDirectory(runDirectory));
    }
}