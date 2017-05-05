package com.hartwig.hmftools.common.context;

import static org.junit.Assert.assertNotNull;

import java.io.File;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.HartwigException;

import org.junit.Test;

public class ProductionRunContextFactoryTest {

    private static final String RESOURCE_DIR = Resources.getResource(
            "context" + File.separator + "ProductionRunContextFactory").getPath();

    @Test
    public void resolveFromRunDirectoryIfNoMetaDataPresent() throws HartwigException {
        final String runDirectory =
                RESOURCE_DIR + File.separator + "170101_HMFregCPCT_FR10002000_FR20003000_CPCT12345678";

        assertNotNull(ProductionRunContextFactory.fromRunDirectory(runDirectory));
    }

    @Test
    public void picksMetadataOverRunNameResolving() throws HartwigException {
        final String runDirectory =
                RESOURCE_DIR + File.separator + "170202_HMFregCPCT_FR10002000_FR20003000_CPCT12345678";

        assertNotNull(ProductionRunContextFactory.fromRunDirectory(runDirectory));
    }
}