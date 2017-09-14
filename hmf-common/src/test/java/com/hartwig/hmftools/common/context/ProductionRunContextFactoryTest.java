package com.hartwig.hmftools.common.context;

import static org.junit.Assert.assertNotNull;

import java.io.File;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.HartwigException;

import org.junit.Test;

public class ProductionRunContextFactoryTest {

    private static final String RESOURCE_DIR = Resources.getResource("context").getPath();

    @Test(expected = HartwigException.class)
    public void throwExceptionWhenNoMetaData() throws HartwigException {
        final String runDirectory = RESOURCE_DIR + File.separator + "DoesNotExist";

        assertNotNull(ProductionRunContextFactory.fromRunDirectory(runDirectory));
    }

    @Test
    public void picksMetadataWhenAvailable() throws HartwigException {
        final String runDirectory = RESOURCE_DIR + File.separator + "RunDirSomatic";

        assertNotNull(ProductionRunContextFactory.fromRunDirectory(runDirectory));
    }
}