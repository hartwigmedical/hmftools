package com.hartwig.hmftools.orange.report;

import java.io.IOException;

import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.OrangeConfigTestFactory;
import com.hartwig.hmftools.orange.OrangeReportTestFactory;
import com.hartwig.hmftools.orange.algo.OrangeAlgo;
import com.hartwig.hmftools.orange.algo.OrangeReport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReportWriterTest {

    @Test
    public void canGenerateTestReportFromTestResources() throws IOException {
        OrangeConfig config = OrangeConfigTestFactory.createDNAConfig();
        OrangeReport report = OrangeAlgo.fromConfig(config).run(config);

        ReportWriter writer = ReportWriterFactory.createInMemoryWriter(config);

        writer.write(report);
    }

    @Test
    public void canGenerateTestReportFromMinimalTestData() throws IOException {
        OrangeReport report = OrangeReportTestFactory.createMinimalTestReport();

        ReportWriter writer = new ReportWriter(false, null, withGermlineReporting());

        writer.write(report);
    }

    @Test
    public void canGenerateTestReportFromProperTestData() throws IOException {
        OrangeReport report = OrangeReportTestFactory.createProperTestReport();

        ReportWriter writer = new ReportWriter(false, null, withGermlineReporting());

        writer.write(report);
    }

    @NotNull
    private static ReportConfig withGermlineReporting() {
        return ImmutableReportConfig.builder().limitJsonOutput(false).reportGermline(true).build();
    }
}
