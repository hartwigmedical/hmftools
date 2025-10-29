package com.hartwig.hmftools.panelbuilder.predesign;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionCenteredAt;

import java.util.List;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.panelbuilder.CustomRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// TODO: about 20% of the existing panel MSI sites are rejected. figure out how to proceed
// Regions around microsatellite instability sites.
public class MsiSites
{
    private static final String FILE_NAME = "msi_regions.tsv";
    private static final int REGION_SIZE = 20;
    private static final String EXTRA_INFO = "MSI";

    private static final Logger LOGGER = LogManager.getLogger(MsiSites.class);

    public static void generateRegions(final String msiSitesFile, final String outputDir)
    {
        LOGGER.info("Generating MSI regions");

        // Don't bother checking overlaps between MSI probes because there are few sites which are far apart.

        List<BasePosition> msiSites = loadMsiSitesFile(msiSitesFile);
        List<CustomRegion> regions = generateRegions(msiSites);
        CustomRegion.writeToFile(regions, outputDir + FILE_NAME);

        LOGGER.info("Done generating MSI regions");
    }

    private static List<BasePosition> loadMsiSitesFile(final String path)
    {
        LOGGER.debug("Loading MSI sites file: {}", path);

        try(DelimFileReader reader = new DelimFileReader(path))
        {
            List<BasePosition> sites = reader.stream().map(row ->
            {
                String chromosome = row.get(FLD_CHROMOSOME);
                int position = row.getInt(FLD_POSITION);
                return new BasePosition(chromosome, position);
            }).toList();

            LOGGER.debug("Loaded {} MSI sites from {}", sites.size(), path);
            return sites;
        }
    }

    private static List<CustomRegion> generateRegions(final List<BasePosition> msiSites)
    {
        return msiSites.stream()
                .map(MsiSites::generateRegion)
                .toList();
    }

    private static CustomRegion generateRegion(final BasePosition msiSite)
    {
        return new CustomRegion(regionCenteredAt(msiSite, REGION_SIZE), EXTRA_INFO);
    }
}
