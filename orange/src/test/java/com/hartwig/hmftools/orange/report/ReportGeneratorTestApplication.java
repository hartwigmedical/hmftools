package com.hartwig.hmftools.orange.report;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.datamodel.cohort.Evaluation;
import com.hartwig.hmftools.datamodel.cohort.ImmutableEvaluation;
import com.hartwig.hmftools.datamodel.isofox.ImmutableIsofoxRecord;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxSvAnnotation;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangeRecord;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.orange.PercentileType;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.orange.ImmutableOrangeConfig;
import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.TestOrangeConfigFactory;
import com.hartwig.hmftools.orange.TestOrangeReportFactory;
import com.hartwig.hmftools.orange.algo.OrangeAlgo;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ReportGeneratorTestApplication {

    private static final Logger LOGGER = LogManager.getLogger(ReportGeneratorTestApplication.class);

    private static final String REPORT_BASE_DIR = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";

    private static final boolean USE_MOCK_DATA_FOR_REPORT = false;
    private static final boolean REMOVE_UNREPORTED_VARIANTS = false;

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
    private static OrangeRecord buildReport(@NotNull OrangeConfig config) throws IOException {
        if (USE_MOCK_DATA_FOR_REPORT) {
            return TestOrangeReportFactory.createProperTestReport();
        }

        OrangeRecord report = OrangeAlgo.fromConfig(config).run(config);

        OrangeRecord withPercentiles = overwriteCohortPercentiles(report);

        OrangeRecord filtered;
        if (REMOVE_UNREPORTED_VARIANTS) {
            filtered = removeUnreported(withPercentiles);
        } else {
            filtered = withPercentiles;
        }

        OrangeRecord finalReport = ImmutableOrangeRecord.builder().from(filtered).sampleId("Test").build();

        return finalReport;
    }

    @NotNull
    private static OrangeRecord overwriteCohortPercentiles(@NotNull OrangeRecord report) {
        // Need to overwrite percentiles since test code doesn't have access to real production cohort percentile files.
        Map<PercentileType, Evaluation> evaluations = Maps.newHashMap();
        evaluations.put(PercentileType.SV_TMB,
                ImmutableEvaluation.builder().cancerType("Skin").panCancerPercentile(0.22).cancerTypePercentile(0.34).build());

        return ImmutableOrangeRecord.builder().from(report).cohortEvaluations(evaluations).build();
    }

    @NotNull
    private static OrangeConfig buildConfig() {
        return ImmutableOrangeConfig.builder()
                .from(TestOrangeConfigFactory.createWGSConfigTumorNormal())
                .outputDir(REPORT_BASE_DIR)
                .addDisclaimer(true)
                .build();
    }

    @NotNull
    private static OrangeRecord removeUnreported(@NotNull OrangeRecord report) {
        ImmutableOrangeRecord.Builder builder = ImmutableOrangeRecord.builder()
                .from(report)
                .purple(ImmutablePurpleRecord.builder()
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
                        .build())
                .linx(ImmutableLinxRecord.builder()
                        .from(report.linx())
                        .allSomaticStructuralVariants(retainReportableStructuralVariants(report.linx().allSomaticStructuralVariants(),
                                report.linx().reportableSomaticBreakends()))
                        .allSomaticFusions(report.linx().reportableSomaticFusions())
                        .additionalSuspectSomaticFusions(Lists.newArrayList())
                        .allSomaticBreakends(report.linx().reportableSomaticBreakends())
                        .additionalSuspectSomaticBreakends(Lists.newArrayList())
                        .allGermlineStructuralVariants(retainReportableStructuralVariants(report.linx().allGermlineStructuralVariants(),
                                report.linx().reportableGermlineBreakends()))
                        .allGermlineBreakends(report.linx().reportableGermlineBreakends())
                        .build());

        if (report.isofox() != null) {
            builder.isofox(ImmutableIsofoxRecord.builder()
                    .from(report.isofox())
                    .allGeneExpressions(Lists.newArrayList())
                    .allFusions(Lists.newArrayList())
                    .allNovelSpliceJunctions(Lists.newArrayList())
                    .build());
        }

        return builder.build();
    }

    @NotNull
    private static List<PurpleGeneCopyNumber> retainReportableCopyNumbers(@NotNull List<PurpleGeneCopyNumber> geneCopyNumbers,
            @NotNull List<PurpleDriver> drivers) {
        List<String> copyNumberDriverGenes = Lists.newArrayList();
        for (PurpleDriver driver : drivers) {
            if (driver.driver() == PurpleDriverType.AMP || driver.driver() == PurpleDriverType.PARTIAL_AMP || driver.driver() == PurpleDriverType.DEL
                    || driver.driver() == PurpleDriverType.GERMLINE_DELETION) {
                copyNumberDriverGenes.add(driver.gene());
            }
        }

        List<PurpleGeneCopyNumber> reportable = Lists.newArrayList();
        for (PurpleGeneCopyNumber geneCopyNumber : geneCopyNumbers) {
            if (copyNumberDriverGenes.contains(geneCopyNumber.geneName())) {
                reportable.add(geneCopyNumber);
            }
        }
        return reportable;
    }

    @Nullable
    private static List<LinxSvAnnotation> retainReportableStructuralVariants(@Nullable List<LinxSvAnnotation> structuralVariants,
            @Nullable List<LinxBreakend> reportableBreakends) {
        if (structuralVariants == null || reportableBreakends == null) {
            return null;
        }

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
