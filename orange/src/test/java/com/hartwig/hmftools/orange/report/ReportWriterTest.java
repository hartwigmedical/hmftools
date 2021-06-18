package com.hartwig.hmftools.orange.report;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.OrangeTestFactory;
import com.hartwig.hmftools.orange.algo.OrangeAlgo;
import com.hartwig.hmftools.orange.algo.OrangeReport;

import org.junit.Test;

public class ReportWriterTest {

    private static final boolean WRITE_TO_PDF = false;

    private static final String REPORT_BASE_DIR = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";

    @Test
    public void canWriteTestReport() throws IOException {
        // TODO Make real test data
        OrangeConfig config = OrangeTestFactory.createTestOrangeConfig();
        OrangeReport report = OrangeAlgo.fromConfig(config).run(config);

        ReportWriter writer = new ReportWriter(REPORT_BASE_DIR);

        writer.write(report, WRITE_TO_PDF);
    }
}