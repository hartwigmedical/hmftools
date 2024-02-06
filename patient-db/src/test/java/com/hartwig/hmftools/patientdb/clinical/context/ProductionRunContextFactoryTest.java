package com.hartwig.hmftools.patientdb.clinical.context;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class ProductionRunContextFactoryTest
{

    private static final String RESOURCE_DIR = Resources.getResource("context").getPath();

    @Test
    public void picksMetadataWhenAvailable() throws IOException
    {
        String runDirectory = RESOURCE_DIR + File.separator + "RunDirSomatic";
        assertNotNull(ProductionRunContextFactory.fromRunDirectory(runDirectory, ""));
    }

    @Test(expected = IOException.class)
    public void throwExceptionWhenNoMetaData() throws IOException
    {
        String runDirectory = RESOURCE_DIR + File.separator + "DoesNotExist";
        assertNotNull(ProductionRunContextFactory.fromRunDirectory(runDirectory, ""));
    }
}