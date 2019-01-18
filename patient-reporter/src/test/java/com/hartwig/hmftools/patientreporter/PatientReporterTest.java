package com.hartwig.hmftools.patientreporter;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testBaseReportData;
import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testSequencedReportData;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.patientreporter.loadStructuralVariants.SvAnalyzerModel;

import org.junit.Test;

public class PatientReporterTest {

    private static final String RUN_DIRECTORY = Resources.getResource("test_run").getPath();

    private static final String FUSION_FILE =
            Resources.getResource("loadStructuralVariants/svAnalysis/CPCT11111111T_fusions.csv").getPath();
    private static final String DISRUPTION_FILE =
            Resources.getResource("loadStructuralVariants/svAnalysis/CPCT11111111T_disruptions.csv").getPath();

    @Test
    public void canRunOnRunDirectory() throws IOException {
        final BaseReportData baseReportData = testBaseReportData();
        final SequencedReportData reporterData = testSequencedReportData();
        final SvAnalyzerModel svAnalyzerModel = SvAnalyzerModel.fromFiles(FUSION_FILE, DISRUPTION_FILE);
        final PatientReporter reporter = ImmutablePatientReporter.of(baseReportData, reporterData, svAnalyzerModel);

        assertNotNull(reporter.run(RUN_DIRECTORY, true, null));
    }
}
