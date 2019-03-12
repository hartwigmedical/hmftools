package com.hartwig.hmftools.patientreporter.cfreport;

import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.ExampleAnalysisTestFactory;

import org.junit.Test;

public class CFReportWriterTest {

    @Test
    public void canGeneratePatientReportForCOLO829() {
        AnalysedPatientReport colo829Report = ExampleAnalysisTestFactory.buildCOLO829();

        CFReportWriter reportWriter = new CFReportWriter();

        reportWriter.writeAnalysedPatientReport(colo829Report, "this_path_does_not_exist");
    }
}
