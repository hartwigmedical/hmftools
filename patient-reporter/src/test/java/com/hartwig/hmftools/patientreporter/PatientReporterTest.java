package com.hartwig.hmftools.patientreporter;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testBaseReportData;
import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testSequencedReportData;
import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testSvAnalyzerModel;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.patientreporter.loadStructuralVariants.SvAnalyzerModel;

import org.junit.Test;

public class PatientReporterTest {

    private static final String RUN_DIRECTORY = Resources.getResource("test_run").getPath();

    @Test
    public void canRunOnRunDirectory() throws IOException {
        final BaseReportData baseReportData = testBaseReportData();
        final SequencedReportData reporterData = testSequencedReportData();
        final SvAnalyzerModel svAnalyzerModel = testSvAnalyzerModel();
        final PatientReporter reporter = ImmutablePatientReporter.of(baseReportData, reporterData, svAnalyzerModel);

        assertNotNull(reporter.run(RUN_DIRECTORY, true, null));
    }
}
