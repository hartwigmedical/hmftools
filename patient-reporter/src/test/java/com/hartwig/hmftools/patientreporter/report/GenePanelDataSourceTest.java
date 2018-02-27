package com.hartwig.hmftools.patientreporter.report;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testHmfReporterData;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.patientreporter.HmfReporterData;
import com.hartwig.hmftools.patientreporter.report.data.GenePanelDataSource;

import org.junit.Test;

import net.sf.jasperreports.engine.JRDataSource;

public class GenePanelDataSourceTest {

    @Test
    public void canCreateGenePanelFor3Genes() throws IOException, HartwigException {
        final HmfReporterData reporterData = testHmfReporterData();
        final JRDataSource dataSource = GenePanelDataSource.fromHmfReporterData(reporterData);
        assertNotNull(dataSource);
    }
}
