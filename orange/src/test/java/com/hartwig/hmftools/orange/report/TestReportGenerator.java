package com.hartwig.hmftools.orange.report;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.linx.ImmutableLinxData;
import com.hartwig.hmftools.common.purple.ImmutablePurpleData;
import com.hartwig.hmftools.orange.ImmutableOrangeConfig;
import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.OrangeConfigTestFactory;
import com.hartwig.hmftools.orange.OrangeReportTestFactory;
import com.hartwig.hmftools.orange.algo.ImmutableOrangeReport;
import com.hartwig.hmftools.orange.algo.OrangeAlgo;
import com.hartwig.hmftools.orange.algo.OrangeReport;
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

    private static final boolean USE_MOCK_DATA_FOR_REPORT = true;

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

        OrangeReport withoutReported = removeUnreported(withPercentiles);

        OrangeReport finalReport = ImmutableOrangeReport.builder().from(withoutReported).sampleId("Test").build();

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
        return ImmutableOrangeConfig.builder().from(OrangeConfigTestFactory.createTestOrangeConfig()).outputDir(REPORT_BASE_DIR).build();
    }

    @NotNull
    private static OrangeReport removeUnreported(@NotNull OrangeReport report) {
        return ImmutableOrangeReport.builder()
                .from(report)
                .purple(ImmutablePurpleData.builder()
                        .from(report.purple())
                        .unreportedGermlineVariants(Lists.newArrayList())
                        .unreportedGainsLosses(Lists.newArrayList())
                        .unreportedSomaticVariants(Lists.newArrayList())
                        .build())
                .linx(ImmutableLinxData.builder().from(report.linx()).unreportedFusions(Lists.newArrayList()).build())
                .build();
    }
}
