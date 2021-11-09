package com.hartwig.hmftools.orange.report;

import java.io.File;
import java.io.IOException;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.linx.ImmutableLinxData;
import com.hartwig.hmftools.common.purple.ImmutablePurpleData;
import com.hartwig.hmftools.orange.ImmutableOrangeConfig;
import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.OrangeTestFactory;
import com.hartwig.hmftools.orange.algo.ImmutableOrangeReport;
import com.hartwig.hmftools.orange.algo.OrangeAlgo;
import com.hartwig.hmftools.orange.algo.OrangeReport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class TestReportGenerator {

    private static final Logger LOGGER = LogManager.getLogger(TestReportGenerator.class);

    private static final String REPORT_BASE_DIR = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";

    public static void main(String[] args) throws IOException {
        OrangeConfig config =
                ImmutableOrangeConfig.builder().from(OrangeTestFactory.createTestOrangeConfig()).outputDir(REPORT_BASE_DIR).build();
        OrangeReport report = OrangeAlgo.fromConfig(config).run(config);

        ReportWriter writer = ReportWriterFactory.createToDiskWriter(config);

        OrangeReport reportWithTestSampleId = removeUnreported(ImmutableOrangeReport.builder().from(report).sampleId("Test").build());

        if (!new File(REPORT_BASE_DIR).isDirectory()) {
            LOGGER.warn("{} is not a directory. Can't write to disk", REPORT_BASE_DIR);
        } else {
            writer.write(reportWithTestSampleId);
        }
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
