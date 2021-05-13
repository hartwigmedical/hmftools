package com.hartwig.hmftools.patientreporter.algo;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.patientreporter.PatientReporterConfig;
import com.hartwig.hmftools.patientreporter.PatientReporterTestFactory;

import org.junit.Test;

public class GenomicAnalyzerTest {

    @Test
    public void canRunOnTestRun() throws IOException {
        AnalysedReportData testReportData = PatientReporterTestFactory.loadTestAnalysedReportData();

        GenomicAnalyzer analyzer = new GenomicAnalyzer(testReportData.germlineReportingModel(),
                testReportData.virusDbModel(),
                testReportData.virusSummaryModel(),
                testReportData.virusBlackListModel());

        PatientReporterConfig config = PatientReporterTestFactory.createTestReporterConfig();

        assertNotNull(analyzer.run("sample", config, LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION));
    }
}
