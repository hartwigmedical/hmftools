package com.hartwig.hmftools.common.linx;

import java.io.IOException;
import java.util.List;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class LinxDataLoader {

    private LinxDataLoader()
    {
    }

    @NotNull
    public static LinxData load(final String tumorSample, final String linxSomaticDir, @Nullable String linxGermlineDir)
            throws IOException
    {
        final String somaticSvAnnotationFile = LinxSvAnnotation.generateFilename(linxSomaticDir, tumorSample, false);
        final String somaticBreakendFile = LinxBreakend.generateFilename(linxSomaticDir, tumorSample);
        final String somaticFusionFile = LinxFusion.generateFilename(linxSomaticDir, tumorSample);
        final String somaticDriversFile = LinxDriver.generateFilename(linxSomaticDir, tumorSample);
        final String somaticDriverCatalogFile = LinxDriver.generateCatalogFilename(linxSomaticDir, tumorSample, true);

        String germlineSvAnnotationFile = null;
        String germlineBreakendFile = null;
        String germlineDisruptionFile = null;
        String germlineDriverCatalogFile = null;

        if (linxGermlineDir != null)
        {
            germlineSvAnnotationFile = LinxSvAnnotation.generateFilename(linxGermlineDir, tumorSample, true);
            germlineBreakendFile = LinxBreakend.generateFilename(linxGermlineDir, tumorSample, true);
            germlineDisruptionFile = LinxGermlineSv.generateFilename(linxGermlineDir, tumorSample);
            germlineDriverCatalogFile = LinxDriver.generateCatalogFilename(linxGermlineDir, tumorSample, false);
        }

        return load(somaticSvAnnotationFile,
                somaticFusionFile,
                somaticBreakendFile,
                somaticDriverCatalogFile,
                somaticDriversFile,
                germlineSvAnnotationFile,
                germlineBreakendFile,
                germlineDisruptionFile,
                germlineDriverCatalogFile);
    }

    @NotNull
    private static LinxData load(@NotNull String somaticStructuralVariantTsv, @NotNull String somaticFusionTsv, @NotNull String somaticBreakendTsv,
            @NotNull String somaticDriverCatalogTsv, @NotNull String somaticDriverTsv, @Nullable String germlineStructuralVariantTsv,
            @Nullable String germlineBreakendTsv, @Nullable String germlineDisruptionTsv, @Nullable String germlineDriverCatalogTsv)
            throws IOException
    {
        List<LinxSvAnnotation> allSomaticStructuralVariants = LinxSvAnnotation.read(somaticStructuralVariantTsv);
        List<LinxDriver> allSomaticDrivers = LinxDriver.read(somaticDriverTsv);

        List<LinxFusion> allSomaticFusions = LinxFusion.read(somaticFusionTsv);
        List<LinxFusion> reportableSomaticFusions = selectReportableFusions(allSomaticFusions);

        List<LinxBreakend> allSomaticBreakends = LinxBreakend.read(somaticBreakendTsv);
        List<LinxBreakend> reportableSomaticBreakends = selectReportableBreakends(allSomaticBreakends);

        List<HomozygousDisruption> somaticHomozygousDisruptions =
                HomozygousDisruptionFactory.extractSomaticFromLinxDriverCatalogTsv(somaticDriverCatalogTsv);

        List<LinxSvAnnotation> allGermlineStructuralVariants = null;
        if (germlineStructuralVariantTsv != null)
        {
            allGermlineStructuralVariants = LinxSvAnnotation.read(germlineStructuralVariantTsv);
        }

        List<LinxBreakend> allGermlineBreakends = null;
        List<LinxBreakend> reportableGermlineBreakends = null;
        if (germlineBreakendTsv != null)
        {
            allGermlineBreakends = LinxBreakend.read(germlineBreakendTsv);
            reportableGermlineBreakends = selectReportableBreakends(allGermlineBreakends);
        }

        List<LinxGermlineSv> allGermlineDisruptions = null;
        List<LinxGermlineSv> reportableGermlineDisruptions = null;
        if (germlineDisruptionTsv != null)
        {
            allGermlineDisruptions = LinxGermlineSv.read(germlineDisruptionTsv);
            reportableGermlineDisruptions = selectReportableGermlineSvs(allGermlineDisruptions, reportableGermlineBreakends);
        }

        List<HomozygousDisruption> germlineHomozygousDisruptions = null;
        if (germlineDriverCatalogTsv != null)
        {
            germlineHomozygousDisruptions = HomozygousDisruptionFactory.extractGermlineFromLinxDriverCatalogTsv(germlineDriverCatalogTsv);
        }

        return ImmutableLinxData.builder()
                .allSomaticStructuralVariants(allSomaticStructuralVariants)
                .somaticDrivers(allSomaticDrivers)
                .allSomaticFusions(allSomaticFusions)
                .reportableSomaticFusions(reportableSomaticFusions)
                .allSomaticBreakends(allSomaticBreakends)
                .reportableSomaticBreakends(reportableSomaticBreakends)
                .somaticHomozygousDisruptions(somaticHomozygousDisruptions)
                .allGermlineStructuralVariants(allGermlineStructuralVariants)
                .allGermlineBreakends(allGermlineBreakends)
                .reportableGermlineBreakends(reportableGermlineBreakends)
                .allGermlineDisruptions(allGermlineDisruptions)
                .reportableGermlineDisruptions(reportableGermlineDisruptions)
                .germlineHomozygousDisruptions(germlineHomozygousDisruptions)
                .build();
    }

    @NotNull
    private static List<LinxFusion> selectReportableFusions(@NotNull List<LinxFusion> fusions)
    {
        List<LinxFusion> reportableFusions = Lists.newArrayList();
        for (LinxFusion fusion : fusions)
        {
            if (fusion.reported())
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
        for (LinxBreakend breakend : breakends)
        {
            if (breakend.reportedDisruption())
            {
                reportableBreakends.add(breakend);
            }
        }
        return reportableBreakends;
    }

    @NotNull
    private static List<LinxGermlineSv> selectReportableGermlineSvs(
            @NotNull List<LinxGermlineSv> germlineSvs, @NotNull List<LinxBreakend> reportableGermlineBreakends) {
        List<LinxGermlineSv> reportableGermlineSvs = Lists.newArrayList();

        if(reportableGermlineBreakends == null)
            return reportableGermlineSvs;

        for (LinxGermlineSv germlineSv : germlineSvs) {

            if(reportableGermlineBreakends.stream().anyMatch(x -> x.svId() == germlineSv.SvId))
                reportableGermlineSvs.add(germlineSv);
        }
        return reportableGermlineSvs;
    }
}
