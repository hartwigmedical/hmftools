package com.hartwig.hmftools.orange.algo.linx;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxHomozygousDisruption;
import com.hartwig.hmftools.orange.conversion.LinxConversion;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class HomozygousDisruptionFactoryTest
{
    private static final String SOMATIC_DRIVERS_CATALOG_TSV = Resources.getResource("linx/sample.linx.driver.catalog.tsv").getPath();
    private static final String GERMLINE_DRIVERS_CATALOG_TSV =
            Resources.getResource("linx/sample.linx.germline.driver.catalog.tsv").getPath();

    /*
    @Test
    public void canExtractSomaticHomozygousDisruptions() throws IOException
    {
        List<HomozygousDisruption> homozygousDisruptions = HomozygousDisruptionFactory.extractSomaticFromLinxDriverCatalogTsv(SOMATIC_DRIVERS_CATALOG_TSV);

        assertEquals(1, homozygousDisruptions.size());

        HomozygousDisruption homozygousDisruption1 = homozygousDisruptions.get(0);
        assertEquals("9", homozygousDisruption1.chromosome());
        assertEquals("p23-p24.1", homozygousDisruption1.chromosomeBand());
        assertEquals("PTPRD", homozygousDisruption1.gene());
    }
    */

    @NotNull
    public static ImmutableHomozygousDisruption.Builder homozygousDisruptionBuilder()
    {
        return ImmutableHomozygousDisruption.builder()
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .gene(Strings.EMPTY)
                .transcript(Strings.EMPTY)
                .isCanonical(false);
    }

    @NotNull
    public static ImmutableLinxHomozygousDisruption.Builder linxHomozygousDisruptionBuilder()
    {
        return ImmutableLinxHomozygousDisruption.builder()
                .from(LinxConversion.convert(homozygousDisruptionBuilder().build()));
    }
}