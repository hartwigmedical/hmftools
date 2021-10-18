package com.hartwig.hmftools.orange.report;

import java.io.File;
import java.io.IOException;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.linx.ImmutableLinxData;
import com.hartwig.hmftools.common.purple.ImmutablePurpleData;
import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.OrangeTestFactory;
import com.hartwig.hmftools.orange.algo.ImmutableOrangeReport;
import com.hartwig.hmftools.orange.algo.OrangeAlgo;
import com.hartwig.hmftools.orange.algo.OrangeReport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReportWriterTest {

    private static final Logger LOGGER = LogManager.getLogger(ReportWriter.class);

    private static final boolean WRITE_TO_DISK = false;
    private static final boolean REMOVE_UNREPORTED = false;

    private static final String REPORT_BASE_DIR = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";

    @Test
    public void canWriteTestReport() throws IOException {
        OrangeConfig config = OrangeTestFactory.createTestOrangeConfig();
        OrangeReport report = OrangeAlgo.fromConfig(config).run(config);

        ReportWriter writer = new ReportWriter(WRITE_TO_DISK, REPORT_BASE_DIR, config.reportConfig());


        OrangeReport reportWithTestSampleId = ImmutableOrangeReport.builder().from(report).sampleId("Test").build();
        OrangeReport finalReport = REMOVE_UNREPORTED ? removeUnreported(reportWithTestSampleId) : reportWithTestSampleId;

        if (WRITE_TO_DISK) {
            if (!new File(REPORT_BASE_DIR).isDirectory()) {
                LOGGER.warn("{} is not a directory. Can't write to disk", REPORT_BASE_DIR);
            } else {
                writer.write(finalReport);
            }
        } else {
            writer.write(finalReport);
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
