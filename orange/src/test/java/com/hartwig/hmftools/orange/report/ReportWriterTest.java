package com.hartwig.hmftools.orange.report;

import java.io.IOException;

import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.TestOrangeConfigFactory;
import com.hartwig.hmftools.orange.TestOrangeReportFactory;
import com.hartwig.hmftools.orange.algo.OrangeAlgo;
import com.hartwig.hmftools.orange.algo.OrangeReport;

import org.junit.Test;

public class ReportWriterTest {

    @Test
    public void canGenerateTestReportForTumorNormalFromTestResources() throws IOException {
        OrangeConfig config = TestOrangeConfigFactory.createDNAConfigTumorNormal();
        OrangeReport report = OrangeAlgo.fromConfig(config).run(config);

        ReportWriter writer = ReportWriterFactory.createInMemoryWriter();

        writer.write(report);
    }

    @Test
    public void canGenerateTestReportForTumorOnlyFromTestResources() throws IOException {
        OrangeConfig config = TestOrangeConfigFactory.createDNAConfigTumorOnly();
        OrangeReport report = OrangeAlgo.fromConfig(config).run(config);

        ReportWriter writer = ReportWriterFactory.createInMemoryWriter();

        writer.write(report);
    }

    @Test
    public void canGenerateTestReportFromMinimalTestData() throws IOException {
        OrangeReport report = TestOrangeReportFactory.createMinimalTestReport();

        ReportWriter writer = ReportWriterFactory.createInMemoryWriter();

        writer.write(report);
    }

    @Test
    public void canGenerateTestReportFromProperTestData() throws IOException {
        OrangeReport report = TestOrangeReportFactory.createProperTestReport();

        ReportWriter writer = ReportWriterFactory.createInMemoryWriter();

        writer.write(report);
    }
}
