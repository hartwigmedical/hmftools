package com.hartwig.hmftools.common.linx;

import java.io.IOException;
import java.util.List;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class LinxDataLoader {
    private LinxDataLoader() {
    }

    @NotNull
    public static LinxData load(final String tumorSample, final String linxSomaticDir, @Nullable String linxGermlineDir)
            throws IOException {
        final String somaticSvAnnotationFile = LinxSvAnnotation.generateFilename(linxSomaticDir, tumorSample, false);
        final String somaticBreakendFile = LinxBreakend.generateFilename(linxSomaticDir, tumorSample);
        final String somaticFusionFile = LinxFusion.generateFilename(linxSomaticDir, tumorSample);
        final String somaticDriversFile = LinxDriver.generateFilename(linxSomaticDir, tumorSample);
        final String somaticDriverCatalogFile = LinxDriver.generateCatalogFilename(linxSomaticDir, tumorSample, true);

        final String germlineSvAnnotationFile =
                linxGermlineDir != null ? LinxSvAnnotation.generateFilename(linxGermlineDir, tumorSample, true) : null;
        final String germlineBreakendFile =
                linxGermlineDir != null ? LinxBreakend.generateFilename(linxGermlineDir, tumorSample, true) : null;
        final String germlineDisruptionFile =
                linxGermlineDir != null ? LinxGermlineSv.generateFilename(linxGermlineDir, tumorSample) : null;

        return load(somaticSvAnnotationFile,
                somaticFusionFile,
                somaticBreakendFile,
                somaticDriverCatalogFile,
                somaticDriversFile,
                germlineSvAnnotationFile,
                germlineBreakendFile,
                germlineDisruptionFile);
    }

    @NotNull
    private static LinxData load(@NotNull String somaticStructuralVariantTsv, @NotNull String somaticFusionTsv, @NotNull String somaticBreakendTsv,
            @NotNull String somaticDriverCatalogTsv, @NotNull String somaticDriverTsv, @Nullable String germlineStructuralVariantTsv,
            @Nullable String germlineBreakendTsv, @Nullable String germlineDisruptionTsv)
            throws IOException {
        List<LinxSvAnnotation> allSomaticStructuralVariants = LinxSvAnnotation.read(somaticStructuralVariantTsv);
        List<LinxDriver> allSomaticDrivers = LinxDriver.read(somaticDriverTsv);

        List<LinxFusion> allSomaticFusions = LinxFusion.read(somaticFusionTsv);
        List<LinxFusion> reportableSomaticFusions = selectReportableFusions(allSomaticFusions);

        List<LinxBreakend> allSomaticBreakends = LinxBreakend.read(somaticBreakendTsv);
        List<LinxBreakend> reportableSomaticBreakends = selectReportableBreakends(allSomaticBreakends);

        List<HomozygousDisruption> somaticHomozygousDisruptions =
                HomozygousDisruptionFactory.extractFromLinxDriverCatalogTsv(somaticDriverCatalogTsv);

        List<LinxSvAnnotation> allGermlineStructuralVariants = null;
        if (germlineStructuralVariantTsv != null) {
            allGermlineStructuralVariants = LinxSvAnnotation.read(germlineStructuralVariantTsv);
        }

        List<LinxBreakend> allGermlineBreakends = null;
        List<LinxBreakend> reportableGermlineBreakends = null;
        if (germlineBreakendTsv != null) {
            allGermlineBreakends = LinxBreakend.read(germlineBreakendTsv);
            reportableGermlineBreakends = selectReportableBreakends(allGermlineBreakends);
        }

        List<LinxGermlineSv> allGermlineDisruptions = null;
        List<LinxGermlineSv> reportableGermlineDisruptions = null;
        if (germlineDisruptionTsv != null) {
            allGermlineDisruptions = LinxGermlineSv.read(germlineDisruptionTsv);
            reportableGermlineDisruptions = selectReportableGermlineSvs(allGermlineDisruptions);
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
                .build();
    }

    @NotNull
    private static List<LinxFusion> selectReportableFusions(@NotNull List<LinxFusion> fusions) {
        List<LinxFusion> reportableFusions = Lists.newArrayList();
        for (LinxFusion fusion : fusions) {
            if (fusion.reported()) {
                reportableFusions.add(fusion);
            }
        }
        return reportableFusions;
    }

    @NotNull
    private static List<LinxBreakend> selectReportableBreakends(@NotNull List<LinxBreakend> breakends) {
        List<LinxBreakend> reportableBreakends = Lists.newArrayList();
        for (LinxBreakend breakend : breakends) {
            if (breakend.reportedDisruption()) {
                reportableBreakends.add(breakend);
            }
        }
        return reportableBreakends;
    }

    @NotNull
    private static List<LinxGermlineSv> selectReportableGermlineSvs(@NotNull List<LinxGermlineSv> germlineSvs) {
        List<LinxGermlineSv> reportableGermlineSvs = Lists.newArrayList();
        for (LinxGermlineSv germlineSv : germlineSvs) {
            if (germlineSv.Reported) {
                reportableGermlineSvs.add(germlineSv);
            }
        }
        return reportableGermlineSvs;
    }
}
