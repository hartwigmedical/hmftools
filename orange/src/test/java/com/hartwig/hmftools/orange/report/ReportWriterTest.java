package com.hartwig.hmftools.orange.report;

import java.io.IOException;

import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.OrangeConfigTestFactory;
import com.hartwig.hmftools.orange.OrangeReportTestFactory;
import com.hartwig.hmftools.orange.algo.OrangeAlgo;
import com.hartwig.hmftools.orange.algo.OrangeReport;

import org.junit.Test;

public class ReportWriterTest {

    @Test
    public void canGenerateTestReportFromTestResources() throws IOException {
        OrangeConfig config = OrangeConfigTestFactory.createTestOrangeConfig();
        OrangeReport report = OrangeAlgo.fromConfig(config).run(config);

        ReportWriter writer = ReportWriterFactory.createInMemoryWriter(config);

        writer.write(report);
    }

    @Test
    public void canGenerateTestReportFromTestData() throws IOException {
        OrangeReport report = OrangeReportTestFactory.createTestReport();

        ReportWriter writer = new ReportWriter(false, null, ImmutableReportConfig.builder().reportGermline(true).build());

        writer.write(report);
    }
}
