package com.hartwig.hmftools.common.linx;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCatalogFile;
import com.hartwig.hmftools.common.driver.DriverType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class HomozygousDisruptionFactory
{
    private static final Logger LOGGER = LogManager.getLogger(HomozygousDisruptionFactory.class);

    private HomozygousDisruptionFactory()
    {
    }

    @NotNull
    public static List<HomozygousDisruption> extractSomaticFromLinxDriverCatalogTsv(@NotNull String linxDriverCatalogTsv)
            throws IOException
    {
        List<DriverCatalog> linxDriversCatalog = DriverCatalogFile.read(linxDriverCatalogTsv);
        LOGGER.debug(" Loaded {} linx driver catalog records from {}", linxDriversCatalog.size(), linxDriverCatalogTsv);

        List<HomozygousDisruption> homozygousDisruptions = extractSomaticHomozygousDisruptions(linxDriversCatalog);
        LOGGER.debug("  Extracted {} somatic homozygous disruptions from linx drivers", homozygousDisruptions.size());
        return homozygousDisruptions;
    }

    @NotNull
    public static List<HomozygousDisruption> extractGermlineFromLinxDriverCatalogTsv(@NotNull String linxDriverCatalogTsv)
            throws IOException
    {
        List<DriverCatalog> linxDriversCatalog = DriverCatalogFile.read(linxDriverCatalogTsv);
        LOGGER.debug(" Loaded {} linx driver catalog records from {}", linxDriversCatalog.size(), linxDriverCatalogTsv);

        List<HomozygousDisruption> homozygousDisruptions = extractGermlineHomozygousDisruptions(linxDriversCatalog);
        LOGGER.debug("  Extracted {} germline homozygous disruptions from linx drivers", homozygousDisruptions.size());
        return homozygousDisruptions;
    }

    private static List<HomozygousDisruption> extractSomaticHomozygousDisruptions(final List<DriverCatalog> driverCatalog)
    {
        List<HomozygousDisruption> homozygousDisruptions = Lists.newArrayList();

        for(DriverCatalog driver : driverCatalog)
        {
            if(driver.driver() == DriverType.HOM_DUP_DISRUPTION || driver.driver() == DriverType.HOM_DEL_DISRUPTION)
            {
                homozygousDisruptions.add(create(driver));
            }
        }

        return homozygousDisruptions;
    }

    private static List<HomozygousDisruption> extractGermlineHomozygousDisruptions(final List<DriverCatalog> driverCatalog)
    {
        List<HomozygousDisruption> homozygousDisruptions = Lists.newArrayList();

        for(DriverCatalog driver : driverCatalog)
        {
            if(driver.driver() == DriverType.GERMLINE_HOM_DUP_DISRUPTION)
            {
                homozygousDisruptions.add(create(driver));
            }
        }

        return homozygousDisruptions;
    }

    private static HomozygousDisruption create(final DriverCatalog driverCatalog)
    {
        return ImmutableHomozygousDisruption.builder()
                .chromosome(driverCatalog.chromosome())
                .chromosomeBand(driverCatalog.chromosomeBand())
                .gene(driverCatalog.gene())
                .transcript(driverCatalog.transcript())
                .isCanonical(driverCatalog.isCanonical())
                .build();
    }
}
