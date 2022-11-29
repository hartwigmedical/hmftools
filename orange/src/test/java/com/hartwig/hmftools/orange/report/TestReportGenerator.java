package com.hartwig.hmftools.orange.report;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.orange.ImmutableOrangeConfig;
import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.OrangeConfigTestFactory;
import com.hartwig.hmftools.orange.OrangeReportTestFactory;
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

public class TestReportGenerator {

    private static final Logger LOGGER = LogManager.getLogger(TestReportGenerator.class);

    private static final String REPORT_BASE_DIR = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";

    private static final boolean USE_MOCK_DATA_FOR_REPORT = false;
    private static final boolean REMOVE_UNREPORTED_VARIANTS = true;

    public static void main(String[] args) throws IOException {
        Configurator.setRootLevel(Level.DEBUG);

        OrangeConfig config = buildConfig();

        ReportWriter writer = ReportWriterFactory.createToDiskWriter(config);

        OrangeReport report = buildReport(config);

        if (!new File(REPORT_BASE_DIR).isDirectory()) {
            LOGGER.warn("{} is not a directory. Can't write to disk", REPORT_BASE_DIR);
        } else {
            writer.write(report);
        }
    }

    @NotNull
    private static OrangeReport buildReport(@NotNull OrangeConfig config) throws IOException {
        if (USE_MOCK_DATA_FOR_REPORT) {
            return OrangeReportTestFactory.createProperTestReport();
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
        return ImmutableOrangeConfig.builder().from(OrangeConfigTestFactory.createDNAConfigTumorNormal()).outputDir(REPORT_BASE_DIR).build();
    }

    @NotNull
    private static OrangeReport removeUnreported(@NotNull OrangeReport report) {
        ImmutableOrangeReport.Builder builder = ImmutableOrangeReport.builder()
                .from(report)
                .purple(ImmutablePurpleInterpretedData.builder()
                        .from(report.purple())
                        .allSomaticVariants(Lists.newArrayList())
                        .additionalSuspectSomaticVariants(Lists.newArrayList())
                        .allGermlineVariants(Lists.newArrayList())
                        .additionalSuspectGermlineVariants(Lists.newArrayList())
                        .allSomaticGeneCopyNumbers(Lists.newArrayList())
                        .suspectGeneCopyNumbersWithLOH(Lists.newArrayList())
                        .allSomaticGainsLosses(Lists.newArrayList())
                        .nearReportableSomaticGains(Lists.newArrayList())
                        .additionalSuspectSomaticGainsLosses(Lists.newArrayList())
                        .allGermlineDeletions(Lists.newArrayList())
                        .build())
                .linx(ImmutableLinxInterpretedData.builder()
                        .from(report.linx())
                        .allStructuralVariants(Lists.newArrayList())
                        .allFusions(Lists.newArrayList())
                        .additionalSuspectFusions(Lists.newArrayList())
                        .allBreakends(Lists.newArrayList())
                        .additionalSuspectBreakends(Lists.newArrayList())
                        .allGermlineDisruptions(Lists.newArrayList())
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
}
