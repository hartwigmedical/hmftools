package com.hartwig.hmftools.orange.report;

import java.io.IOException;

import com.hartwig.hmftools.orange.ImmutableOrangeConfig;
import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.OrangeConfigTestFactory;
import com.hartwig.hmftools.orange.OrangeReportTestFactory;
import com.hartwig.hmftools.orange.algo.OrangeAlgo;
import com.hartwig.hmftools.orange.algo.OrangeReport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReportWriterTest {

    @Test
    public void canGenerateTestReportForTumorNormalFromTestResources() throws IOException {
        OrangeConfig config = OrangeConfigTestFactory.createDNAConfigTumorNormal();
        OrangeReport report = OrangeAlgo.fromConfig(config).run(config);

        ReportWriter writer = ReportWriterFactory.createInMemoryWriter(config);

        writer.write(report);
    }

    @Test
    public void canGenerateTestReportForTumorOnlyFromTestResources() throws IOException {
        OrangeConfig config = OrangeConfigTestFactory.createDNAConfigTumorOnly();
        OrangeReport report = OrangeAlgo.fromConfig(config).run(config);

        ReportWriter writer = ReportWriterFactory.createInMemoryWriter(config);

        writer.write(report);
    }

    @Test
    public void canGenerateTestReportFromMinimalTestData() throws IOException {
        OrangeReport report = OrangeReportTestFactory.createMinimalTestReport();

        ReportWriter writer = ReportWriterFactory.createInMemoryWriter(noGermlineReporting());

        writer.write(report);
    }

    @Test
    public void canGenerateTestReportFromProperTestData() throws IOException {
        OrangeReport report = OrangeReportTestFactory.createProperTestReport();

        ReportWriter writer = ReportWriterFactory.createInMemoryWriter(withGermlineReporting());

        writer.write(report);
    }

    @NotNull
    private static OrangeConfig withGermlineReporting() {
        return withReportingConfig(ImmutableReportConfig.builder().limitJsonOutput(false).reportGermline(true).build());
    }

    @NotNull
    private static OrangeConfig noGermlineReporting() {
        return withReportingConfig(ImmutableReportConfig.builder().limitJsonOutput(false).reportGermline(false).build());
    }

    @NotNull
    private static OrangeConfig withReportingConfig(@NotNull ReportConfig reportConfig) {
        return ImmutableOrangeConfig.builder()
                .from(OrangeConfigTestFactory.createDNAConfigTumorNormal())
                .reportConfig(reportConfig)
                .build();
    }
}
