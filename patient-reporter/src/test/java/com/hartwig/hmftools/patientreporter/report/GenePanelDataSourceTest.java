package com.hartwig.hmftools.patientreporter.report;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testSequencedReportData;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.hartwig.hmftools.patientreporter.SequencedReportData;
import com.hartwig.hmftools.patientreporter.report.data.GenePanelDataSource;

import org.junit.Test;

import net.sf.jasperreports.engine.JRDataSource;

public class GenePanelDataSourceTest {

    @Test
    public void canCreateGenePanelFor3Genes() throws IOException {
        final SequencedReportData reporterData = testSequencedReportData();
        final JRDataSource dataSource = GenePanelDataSource.fromSequencedReportData(reporterData);
        assertNotNull(dataSource);
    }
}
