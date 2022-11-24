package com.hartwig.hmftools.common.linx;

import java.io.IOException;
import java.util.List;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class LinxDataLoader
{
    private LinxDataLoader() {}

    @NotNull
    public static LinxData load(final String tumorSample, final String linxSomaticDir, @Nullable String linxGermlineDir) throws IOException
    {
        final String svAnnotationFile = LinxSvAnnotation.generateFilename(linxSomaticDir, tumorSample, false);
        final String breakendFile = LinxBreakend.generateFilename(linxSomaticDir, tumorSample);
        final String fusionFile = LinxFusion.generateFilename(linxSomaticDir, tumorSample);
        final String driversFile = LinxDriver.generateFilename(linxSomaticDir, tumorSample);
        final String driverCatalogFile = LinxDriver.generateCatalogFilename(linxSomaticDir, tumorSample, true);
        final String germlineSvFile = linxGermlineDir != null ? LinxGermlineSv.generateFilename(linxGermlineDir, tumorSample) : null;

        return load(svAnnotationFile, fusionFile, breakendFile, driverCatalogFile, driversFile, germlineSvFile);
    }

    @NotNull
    private static LinxData load(@NotNull String linxStructuralVariantTsv, @NotNull String linxFusionTsv, @NotNull String linxBreakendTsv,
            @NotNull String linxDriverCatalogTsv, @NotNull String linxDriverTsv, @Nullable String linxGermlineDisruptionTsv)
            throws IOException
    {
        List<LinxFusion> allFusions = LinxFusion.read(linxFusionTsv);
        List<LinxFusion> reportableFusions = selectReportableFusions(allFusions);

        List<LinxBreakend> allBreakends = LinxBreakend.read(linxBreakendTsv);
        List<LinxBreakend> reportableBreakends = selectReportableBreakends(allBreakends);

        List<HomozygousDisruption> homozygousDisruptions =
                HomozygousDisruptionFactory.extractFromLinxDriverCatalogTsv(linxDriverCatalogTsv);

        List<LinxGermlineSv> allGermlineDisruptions = null;
        List<LinxGermlineSv> reportableGermlineDisruptions = null;
        if(linxGermlineDisruptionTsv != null)
        {
            allGermlineDisruptions = LinxGermlineSv.read(linxGermlineDisruptionTsv);
            reportableGermlineDisruptions = selectReportableGermlineSvs(allGermlineDisruptions);
        }

        return ImmutableLinxData.builder()
                .allStructuralVariants(LinxSvAnnotation.read(linxStructuralVariantTsv))
                .drivers(LinxDriver.read(linxDriverTsv))
                .allFusions(allFusions)
                .reportableFusions(reportableFusions)
                .allBreakends(allBreakends)
                .reportableBreakends(reportableBreakends)
                .homozygousDisruptions(homozygousDisruptions)
                .allGermlineDisruptions(allGermlineDisruptions)
                .reportableGermlineDisruptions(reportableGermlineDisruptions)
                .build();
    }

    @NotNull
    private static List<LinxFusion> selectReportableFusions(@NotNull List<LinxFusion> fusions)
    {
        List<LinxFusion> reportableFusions = Lists.newArrayList();
        for(LinxFusion fusion : fusions)
        {
            if(fusion.reported())
            {
                reportableFusions.add(fusion);
            }
        }
        return reportableFusions;
    }

    @NotNull
    private static List<LinxBreakend> selectReportableBreakends(@NotNull List<LinxBreakend> breakends)
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

    @NotNull
    private static List<LinxGermlineSv> selectReportableGermlineSvs(@NotNull List<LinxGermlineSv> germlineSvs)
    {
        List<LinxGermlineSv> reportableGermlineSvs = Lists.newArrayList();
        for(LinxGermlineSv germlineSv : germlineSvs)
        {
            if(germlineSv.Reported)
            {
                reportableGermlineSvs.add(germlineSv);
            }
        }
        return reportableGermlineSvs;
    }
}
