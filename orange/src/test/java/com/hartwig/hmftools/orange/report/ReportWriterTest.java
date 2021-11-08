package com.hartwig.hmftools.orange.report;

import java.io.IOException;

import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.OrangeTestFactory;
import com.hartwig.hmftools.orange.algo.OrangeAlgo;
import com.hartwig.hmftools.orange.algo.OrangeReport;

import org.junit.Test;

public class ReportWriterTest {

    @Test
    public void canGenerateTestReport() throws IOException {
        OrangeConfig config = OrangeTestFactory.createTestOrangeConfig();
        OrangeReport report = OrangeAlgo.fromConfig(config).run(config);

        ReportWriter writer = ReportWriterFactory.createInMemoryWriter(config);

        writer.write(report);
    }
}
