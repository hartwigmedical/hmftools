package com.hartwig.hmftools.patientreporter.cfreport;

import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.ExampleAnalysisTestFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

import static org.junit.Assert.assertNotNull;

public class CFReportWriterTest {

    private static final String REPORT_BASE_DIR = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";

    @Test
    public void canGeneratePatientReportForCOLO829() throws IOException {
        AnalysedPatientReport colo829Report = ExampleAnalysisTestFactory.buildCOLO829();

        CFReportWriter reportWriter = new CFReportWriter();
        reportWriter.writeAnalysedPatientReport(colo829Report,
                getReportFilePath("hmf_test_sequence_report_" + String.valueOf(System.currentTimeMillis()) + ".pdf"));
    }

    @Test
    public void canGeneratePatientReportForCompletelyFilledInReport() throws IOException {
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledIn();

        CFReportWriter reportWriter = new CFReportWriter();
        reportWriter.writeAnalysedPatientReport(patientReport,
                getReportFilePath("hmf_full_test_sequence_report_" + String.valueOf(System.currentTimeMillis()) + ".pdf"));
    }

    @NotNull
    private static String getReportFilePath(@NotNull String filename) {
        return REPORT_BASE_DIR + File.separator + filename;
    }
}
