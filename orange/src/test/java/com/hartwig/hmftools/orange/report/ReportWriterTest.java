package com.hartwig.hmftools.orange.report;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.OrangeTestFactory;
import com.hartwig.hmftools.orange.algo.ImmutableOrangeReport;
import com.hartwig.hmftools.orange.algo.OrangeAlgo;
import com.hartwig.hmftools.orange.algo.OrangeReport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.Test;

public class ReportWriterTest {

    private static final Logger LOGGER = LogManager.getLogger(ReportWriter.class);

    private static final boolean WRITE_TO_PDF = true;

    private static final String REPORT_BASE_DIR = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";

    @Test
    public void canWriteTestReport() throws IOException {
        OrangeConfig config = OrangeTestFactory.createTestOrangeConfig();
        OrangeReport report = OrangeAlgo.fromConfig(config).run(config);

        ReportWriter writer = new ReportWriter(WRITE_TO_PDF, REPORT_BASE_DIR, config.reportConfig());

        OrangeReport reportWithTestSampleId = ImmutableOrangeReport.builder().from(report).sampleId("Test").build();

        if (WRITE_TO_PDF) {
            if (!new File(REPORT_BASE_DIR).isDirectory()) {
                LOGGER.warn("{} is not a directory. Can't write PDF", REPORT_BASE_DIR);
            } else {
                writer.write(reportWithTestSampleId);
            }
        } else {
            writer.write(reportWithTestSampleId);
        }
    }
}
