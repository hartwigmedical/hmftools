package com.hartwig.hmftools.orange.report;

import java.io.IOException;

import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.TestOrangeConfigFactory;
import com.hartwig.hmftools.orange.TestOrangeReportFactory;
import com.hartwig.hmftools.orange.algo.OrangeAlgo;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReportWriterTest
{
    @Test
    public void canGenerateTestReportForTargetedFromTestResources() throws IOException
    {
        run(TestOrangeConfigFactory.createTargetedConfig());
    }

    @Test
    public void canGenerateTestReportForWGSTumorOnlyFromTestResources() throws IOException
    {
        run(TestOrangeConfigFactory.createWGSConfigTumorOnly());
    }

    @Test
    public void canGenerateTestReportForWGSTumorNormalFromTestResources() throws IOException
    {
        run(TestOrangeConfigFactory.createWGSConfigTumorNormal());
    }

    @Test
    public void canGenerateTestReportFromMinimalTestData() throws IOException
    {
        OrangeRecord report = TestOrangeReportFactory.createMinimalTestReport();

        ReportWriter writer = ReportWriterFactory.createInMemoryWriter();

        writer.write(report);
    }

    @Test
    public void canGenerateTestReportFromProperTestData() throws IOException
    {
        OrangeRecord report = TestOrangeReportFactory.createProperTestReport();

        ReportWriter writer = ReportWriterFactory.createInMemoryWriter();

        writer.write(report);
    }

    private static void run(@NotNull OrangeConfig config) throws IOException
    {
        OrangeAlgo algo = OrangeAlgo.fromConfig(config);
        algo.setSuppressGeneWarnings();
        OrangeRecord report = algo.run(config);

        ReportWriter writer = ReportWriterFactory.createInMemoryWriter();

        writer.write(report);
    }
}
