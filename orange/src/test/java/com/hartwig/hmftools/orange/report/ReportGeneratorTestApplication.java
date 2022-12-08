package com.hartwig.hmftools.orange.report;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.orange.ImmutableOrangeConfig;
import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.TestOrangeConfigFactory;
import com.hartwig.hmftools.orange.TestOrangeReportFactory;
import com.hartwig.hmftools.orange.algo.ImmutableOrangeReport;
import com.hartwig.hmftools.orange.algo.OrangeAlgo;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.algo.isofox.ImmutableIsofoxInterpretedData;
import com.hartwig.hmftools.orange.algo.linx.ImmutableLinxInterpretedData;
import com.hartwig.hmftools.orange.algo.purple.ImmutablePurpleInterpretedData;
import com.hartwig.hmftools.orange.cohort.datamodel.Evaluation;
import com.hartwig.hmftools.orange.cohort.datamodel.ImmutableEvaluation;
import com.hartwig.hmftools.orange.cohort.percentile.PercentileType;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class ReportGeneratorTestApplication {

    private static final Logger LOGGER = LogManager.getLogger(ReportGeneratorTestApplication.class);

    private static final String REPORT_BASE_DIR = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";

    private static final boolean USE_MOCK_DATA_FOR_REPORT = false;
    private static final boolean REMOVE_UNREPORTED_VARIANTS = true;

    public static void main(String[] args) throws IOException {
        Configurator.setRootLevel(Level.DEBUG);

        OrangeConfig config = buildConfig();

        ReportWriter writer = ReportWriterFactory.createToDiskWriter(config);

        if (!new File(REPORT_BASE_DIR).isDirectory()) {
            LOGGER.warn("{} is not a directory. Can't write to disk", REPORT_BASE_DIR);
        } else {
            LOGGER.info("Deleting plot dir");
            deleteDir(new File(REPORT_BASE_DIR + File.separator + "plot"));
            writer.write(buildReport(config));
        }
    }

    @NotNull
    private static OrangeReport buildReport(@NotNull OrangeConfig config) throws IOException {
        if (USE_MOCK_DATA_FOR_REPORT) {
            return TestOrangeReportFactory.createProperTestReport();
        }

        OrangeReport report = OrangeAlgo.fromConfig(config).run(config);

        OrangeReport withPercentiles = overwriteCohortPercentiles(report);

        OrangeReport filtered;
        if (REMOVE_UNREPORTED_VARIANTS) {
            filtered = removeUnreported(withPercentiles);
        } else {
            filtered = withPercentiles;
        }

        OrangeReport finalReport = ImmutableOrangeReport.builder().from(filtered).sampleId("Test").build();

        return finalReport;
    }

    @NotNull
    private static OrangeReport overwriteCohortPercentiles(@NotNull OrangeReport report) {
        // Need to overwrite percentiles since test code doesn't have access to real production cohort percentile files.
        Map<PercentileType, Evaluation> evaluations = Maps.newHashMap();
        evaluations.put(PercentileType.SV_TMB,
                ImmutableEvaluation.builder().cancerType("Skin").panCancerPercentile(0.22).cancerTypePercentile(0.34).build());

        return ImmutableOrangeReport.builder().from(report).cohortEvaluations(evaluations).build();
    }

    @NotNull
    private static OrangeConfig buildConfig() {
        return ImmutableOrangeConfig.builder()
                .from(TestOrangeConfigFactory.createDNAConfigTumorNormal())
                .outputDir(REPORT_BASE_DIR)
                .build();
    }

    @NotNull
    private static OrangeReport removeUnreported(@NotNull OrangeReport report) {
        ImmutableOrangeReport.Builder builder = ImmutableOrangeReport.builder()
                .from(report)
                .purple(ImmutablePurpleInterpretedData.builder()
                        .from(report.purple())
                        .allSomaticVariants(report.purple().reportableSomaticVariants())
                        .additionalSuspectSomaticVariants(Lists.newArrayList())
                        .allGermlineVariants(report.purple().reportableGermlineVariants())
                        .additionalSuspectGermlineVariants(Lists.newArrayList())
                        .allSomaticCopyNumbers(Lists.newArrayList())
                        .allSomaticGeneCopyNumbers(retainReportableCopyNumbers(report.purple().allSomaticGeneCopyNumbers(),
                                report.purple().somaticDrivers()))
                        .suspectGeneCopyNumbersWithLOH(Lists.newArrayList())
                        .allSomaticGainsLosses(report.purple().reportableSomaticGainsLosses())
                        .nearReportableSomaticGains(Lists.newArrayList())
                        .additionalSuspectSomaticGainsLosses(Lists.newArrayList())
                        .allGermlineDeletions(report.purple().reportableGermlineDeletions())
                        .build())
                .linx(ImmutableLinxInterpretedData.builder()
                        .from(report.linx())
                        .allStructuralVariants(retainReportableStructuralVariants(report.linx().allStructuralVariants(),
                                report.linx().reportableBreakends()))
                        .allFusions(report.linx().reportableFusions())
                        .additionalSuspectFusions(Lists.newArrayList())
                        .allBreakends(report.linx().reportableBreakends())
                        .additionalSuspectBreakends(Lists.newArrayList())
                        .allGermlineDisruptions(report.linx().reportableGermlineDisruptions())
                        .build());

        if (report.isofox() != null) {
            builder.isofox(ImmutableIsofoxInterpretedData.builder()
                    .from(report.isofox())
                    .allGeneExpressions(Lists.newArrayList())
                    .allFusions(Lists.newArrayList())
                    .allNovelSpliceJunctions(Lists.newArrayList())
                    .build());
        }

        return builder.build();
    }

    @NotNull
    private static List<GeneCopyNumber> retainReportableCopyNumbers(@NotNull List<GeneCopyNumber> geneCopyNumbers,
            @NotNull List<DriverCatalog> drivers) {
        List<String> copyNumberDriverGenes = Lists.newArrayList();
        for (DriverCatalog driver : drivers) {
            if (driver.driver() == DriverType.AMP || driver.driver() == DriverType.PARTIAL_AMP || driver.driver() == DriverType.DEL
                    || driver.driver() == DriverType.GERMLINE_DELETION) {
                copyNumberDriverGenes.add(driver.gene());
            }
        }

        List<GeneCopyNumber> reportable = Lists.newArrayList();
        for (GeneCopyNumber geneCopyNumber : geneCopyNumbers) {
            if (copyNumberDriverGenes.contains(geneCopyNumber.geneName())) {
                reportable.add(geneCopyNumber);
            }
        }
        return reportable;
    }

    @NotNull
    private static List<LinxSvAnnotation> retainReportableStructuralVariants(@NotNull List<LinxSvAnnotation> structuralVariants,
            @NotNull List<LinxBreakend> reportableBreakends) {
        List<LinxSvAnnotation> reportable = Lists.newArrayList();
        for (LinxSvAnnotation structuralVariant : structuralVariants) {
            if (isReportableSv(structuralVariant, reportableBreakends)) {
                reportable.add(structuralVariant);
            }
        }
        return reportable;
    }

    private static boolean isReportableSv(@NotNull LinxSvAnnotation structuralVariant, @NotNull List<LinxBreakend> reportableBreakends) {
        for (LinxBreakend breakend : reportableBreakends) {
            if (breakend.svId() == structuralVariant.svId()) {
                return true;
            }
        }
        return false;
    }

    private static void deleteDir(@NotNull File file) {
        File[] contents = file.listFiles();
        if (contents != null) {
            for (File content : contents) {
                if (!Files.isSymbolicLink(content.toPath())) {
                    deleteDir(content);
                }
            }
        }
        file.delete();
    }
}
