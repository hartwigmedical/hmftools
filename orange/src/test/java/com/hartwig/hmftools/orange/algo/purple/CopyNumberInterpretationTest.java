package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogTestFactory;
import com.hartwig.hmftools.common.drivercatalog.DriverType;

import org.junit.Test;

public class CopyNumberInterpretationTest {

    @Test
    public void canRenderDisplay() {
        assertEquals("full gain", CopyNumberInterpretation.FULL_GAIN.display());
    }

    @Test
    public void canConvertCNADrivers() {
        DriverCatalog amp = DriverCatalogTestFactory.builder().driver(DriverType.AMP).build();
        assertEquals(CopyNumberInterpretation.FULL_GAIN, CopyNumberInterpretation.fromCNADriver(amp));

        DriverCatalog partialAmp = DriverCatalogTestFactory.builder().driver(DriverType.PARTIAL_AMP).build();
        assertEquals(CopyNumberInterpretation.PARTIAL_GAIN, CopyNumberInterpretation.fromCNADriver(partialAmp));

        DriverCatalog fullLoss = DriverCatalogTestFactory.builder().driver(DriverType.DEL).maxCopyNumber(0).build();
        assertEquals(CopyNumberInterpretation.FULL_LOSS, CopyNumberInterpretation.fromCNADriver(fullLoss));

        DriverCatalog partialLoss = DriverCatalogTestFactory.builder().driver(DriverType.DEL).maxCopyNumber(2).build();
        assertEquals(CopyNumberInterpretation.PARTIAL_LOSS, CopyNumberInterpretation.fromCNADriver(partialLoss));
    }

    @Test (expected = IllegalStateException.class)
    public void crashOnNonCNADriver() {
        CopyNumberInterpretation.fromCNADriver(DriverCatalogTestFactory.builder().driver(DriverType.GERMLINE_MUTATION).build());
    }
}