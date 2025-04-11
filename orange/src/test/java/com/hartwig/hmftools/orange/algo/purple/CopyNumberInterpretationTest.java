package com.hartwig.hmftools.orange.algo.purple;

import static com.hartwig.hmftools.orange.algo.purple.CopyNumberInterpretationUtil.display;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCatalogTestFactory;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;

import org.junit.Test;

public class CopyNumberInterpretationTest
{
    @Test
    public void canRenderDisplay()
    {
        assertEquals("full gain", display(CopyNumberInterpretation.FULL_GAIN));
    }

    @Test
    public void canConvertCNADrivers()
    {
        DriverCatalog amp = DriverCatalogTestFactory.builder().driver(DriverType.AMP).build();
        assertEquals(CopyNumberInterpretation.FULL_GAIN, CopyNumberInterpretationUtil.fromCNADriver(amp));

        DriverCatalog partialAmp = DriverCatalogTestFactory.builder().driver(DriverType.PARTIAL_AMP).build();
        assertEquals(CopyNumberInterpretation.PARTIAL_GAIN, CopyNumberInterpretationUtil.fromCNADriver(partialAmp));

        DriverCatalog fullDel = DriverCatalogTestFactory.builder().driver(DriverType.DEL).maxCopyNumber(0).build();
        assertEquals(CopyNumberInterpretation.FULL_DEL, CopyNumberInterpretationUtil.fromCNADriver(fullDel));

        DriverCatalog partialDel = DriverCatalogTestFactory.builder().driver(DriverType.DEL).maxCopyNumber(2).build();
        assertEquals(CopyNumberInterpretation.PARTIAL_DEL, CopyNumberInterpretationUtil.fromCNADriver(partialDel));
    }

    @Test(expected = IllegalStateException.class)
    public void crashOnNonCNADriver()
    {
        CopyNumberInterpretationUtil.fromCNADriver(DriverCatalogTestFactory.builder().driver(DriverType.GERMLINE_MUTATION).build());
    }
}