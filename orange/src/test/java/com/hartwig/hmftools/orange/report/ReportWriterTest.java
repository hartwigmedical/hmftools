package com.hartwig.hmftools.orange.report;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.OrangeTestFactory;
import com.hartwig.hmftools.orange.algo.ImmutableOrangeReport;
import com.hartwig.hmftools.orange.algo.OrangeAlgo;
import com.hartwig.hmftools.orange.algo.OrangeReport;

import org.junit.Test;

public class ReportWriterTest {

    private static final boolean WRITE_TO_PDF = false;

    private static final String REPORT_BASE_DIR = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";

    @Test
    public void canWriteTestReport() throws IOException {
        ReportWriter writer = new ReportWriter(WRITE_TO_PDF, REPORT_BASE_DIR);

        OrangeConfig config = OrangeTestFactory.createTestOrangeConfig();
        OrangeReport report = OrangeAlgo.fromConfig(config).run(config);

        OrangeReport reportWithTestSampleId = ImmutableOrangeReport.builder().from(report).sampleId("Test").build();

        writer.write(reportWithTestSampleId);
    }
}