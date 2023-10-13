package com.hartwig.hmftools.common.linx;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class HomozygousDisruptionFactoryTest
{
    private static final String SOMATIC_DRIVERS_CATALOG_TSV = Resources.getResource("linx/sample.linx.driver.catalog.tsv").getPath();
    private static final String GERMLINE_DRIVERS_CATALOG_TSV =
            Resources.getResource("linx/sample.linx.germline.driver.catalog.tsv").getPath();

    @Test
    public void canExtractSomaticHomozygousDisruptions() throws IOException
    {
        List<HomozygousDisruption> homozygousDisruptions =
                HomozygousDisruptionFactory.extractSomaticFromLinxDriverCatalogTsv(SOMATIC_DRIVERS_CATALOG_TSV);

        assertEquals(1, homozygousDisruptions.size());

        HomozygousDisruption homozygousDisruption1 = homozygousDisruptions.get(0);
        assertEquals("9", homozygousDisruption1.chromosome());
        assertEquals("p23-p24.1", homozygousDisruption1.chromosomeBand());
        assertEquals("PTPRD", homozygousDisruption1.gene());
    }

    @Test
    public void canExtractGermlineHomozygousDisruptions() throws IOException
    {
        List<HomozygousDisruption> homozygousDisruptions =
                HomozygousDisruptionFactory.extractGermlineFromLinxDriverCatalogTsv(GERMLINE_DRIVERS_CATALOG_TSV);

        assertEquals(1, homozygousDisruptions.size());

        HomozygousDisruption homozygousDisruption1 = homozygousDisruptions.get(0);
        assertEquals("13", homozygousDisruption1.chromosome());
        assertEquals("q13.1", homozygousDisruption1.chromosomeBand());
        assertEquals("BRCA2", homozygousDisruption1.gene());
        assertEquals("ENST00000544455", homozygousDisruption1.transcript());
        assertTrue(homozygousDisruption1.isCanonical());
    }
}