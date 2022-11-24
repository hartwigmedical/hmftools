package com.hartwig.hmftools.common.linx;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class LinxDataLoader
{
    private static final Logger LOGGER = LogManager.getLogger(LinxDataLoader.class);

    private LinxDataLoader() {}

    @NotNull
    public static LinxData load(final String tumorSample, final String linxSomaticDir, final String linxGermlineDir) throws IOException
    {
        final String svAnnotationFile = LinxSvAnnotation.generateFilename(linxSomaticDir, tumorSample, false);
        final String breakendFile = LinxBreakend.generateFilename(linxSomaticDir, tumorSample);
        final String fusionFile = LinxFusion.generateFilename(linxSomaticDir, tumorSample);
        final String driversFile = LinxDriver.generateFilename(linxSomaticDir, tumorSample);
        final String driverCatalogFile = LinxDriver.generateCatalogFilename(linxSomaticDir, tumorSample, true);
        final String germlineSvFile = LinxGermlineSv.generateFilename(linxGermlineDir, tumorSample);

        return load(svAnnotationFile, fusionFile, breakendFile, driverCatalogFile, driversFile, germlineSvFile);
    }

    @NotNull
    private static LinxData load(@Nullable String linxStructuralVariantTsv, @NotNull String linxFusionTsv, @NotNull String linxBreakendTsv,
            @NotNull String linxDriverCatalogTsv, @Nullable String linxDriverTsv, @Nullable String linxGermlineDisruptionTsv)
            throws IOException
    {
        LOGGER.info("Loading LINX data from {}", new File(linxFusionTsv).getParent());
        List<LinxFusion> allFusions = LinxFusion.read(linxFusionTsv);

        List<LinxFusion> reportableFusions = Lists.newArrayList();
        for(LinxFusion fusion : allFusions)
        {
            if(fusion.reported())
            {
                reportableFusions.add(fusion);
            }
        }
        LOGGER.info(" Loaded {} fusions (of which {} are reportable) from {}", allFusions.size(), reportableFusions.size(), linxFusionTsv);

        List<LinxSvAnnotation> allStructuralVariants = Lists.newArrayList();
        if(linxStructuralVariantTsv != null)
        {
            allStructuralVariants = LinxSvAnnotation.read(linxStructuralVariantTsv);
            LOGGER.info(" Loaded {} structural variants from {}", allStructuralVariants.size(), linxStructuralVariantTsv);
        }

        List<LinxBreakend> allBreakends = LinxBreakend.read(linxBreakendTsv);
        List<LinxBreakend> reportableBreakends = selectReportable(allBreakends);
        LOGGER.info(" Loaded {} breakends (of which {} are reportable) from {}",
                allBreakends.size(),
                reportableBreakends.size(),
                linxBreakendTsv);

        List<HomozygousDisruption> homozygousDisruptions =
                HomozygousDisruptionFactory.extractFromLinxDriverCatalogTsv(linxDriverCatalogTsv);
        LOGGER.info(" Loaded {} reportable homozygous disruptions from {}", homozygousDisruptions.size(), linxDriverCatalogTsv);

        List<LinxDriver> allDrivers = Lists.newArrayList();
        if(linxDriverTsv != null)
        {
            allDrivers = LinxDriver.read(linxDriverTsv);
            LOGGER.info(" Loaded {} drivers from {}", allDrivers.size(), linxDriverTsv);
        }

        List<LinxGermlineSv> allGermlineDisruptions = Lists.newArrayList();
        List<LinxGermlineSv> reportableGermlineDisruptions = Lists.newArrayList();
        if(linxGermlineDisruptionTsv != null)
        {
            allGermlineDisruptions = LinxGermlineSv.read(linxGermlineDisruptionTsv);

            for(LinxGermlineSv germlineDisruption : allGermlineDisruptions)
            {
                if(germlineDisruption.Reported)
                {
                    reportableGermlineDisruptions.add(germlineDisruption);
                }
            }
            LOGGER.info(" Loaded {} germline disruptions (of which {} are reportable) from {}",
                    allGermlineDisruptions.size(),
                    reportableGermlineDisruptions.size(),
                    linxGermlineDisruptionTsv);
        }

        return ImmutableLinxData.builder()
                .allStructuralVariants(allStructuralVariants)
                .allFusions(allFusions)
                .reportableFusions(reportableFusions)
                .allBreakends(allBreakends)
                .reportableBreakends(reportableBreakends)
                .homozygousDisruptions(homozygousDisruptions)
                .drivers(allDrivers)
                .allGermlineDisruptions(allGermlineDisruptions)
                .reportableGermlineDisruptions(reportableGermlineDisruptions)
                .build();
    }

    @NotNull
    private static List<LinxBreakend> selectReportable(@NotNull List<LinxBreakend> breakends)
    {
        List<LinxBreakend> reportableBreakends = Lists.newArrayList();
        for(LinxBreakend breakend : breakends)
        {
            if(breakend.reportedDisruption())
            {
                reportableBreakends.add(breakend);
            }
        }
        return reportableBreakends;
    }
}
