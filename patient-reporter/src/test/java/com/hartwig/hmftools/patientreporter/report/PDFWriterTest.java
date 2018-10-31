package com.hartwig.hmftools.patientreporter.report;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testBaseReportData;
import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testSampleReport;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Optional;

import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.ImmutableNotAnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.NotAnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.qcfail.NotAnalysableReason;
import com.hartwig.hmftools.patientreporter.qcfail.NotAnalysableStudy;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.exception.DRException;

public class PDFWriterTest {

    private static final boolean WRITE_TO_PDF = false;

    private static final String REPORT_BASE_DIR = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";

    @Test
    public void canGenerateSequenceReport() throws DRException, IOException {
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.buildCOLO829();

        JasperReportBuilder mainReport = PDFWriter.generatePatientReport(patientReport);
        assertNotNull(mainReport);

        if (WRITE_TO_PDF) {
            mainReport.toPdf(new FileOutputStream(REPORT_BASE_DIR + File.separator + "hmf_test_sequence_report.pdf"));
        }
    }

    @Test
    public void canGenerateLowTumorPercentageReport() throws DRException, IOException {
        JasperReportBuilder report = generateNotAnalysableCPCTReport(0.1, NotAnalysableReason.LOW_TUMOR_PERCENTAGE);
        assertNotNull(report);

        if (WRITE_TO_PDF) {
            report.toPdf(new FileOutputStream(REPORT_BASE_DIR + File.separator + "hmf_low_tumor_percentage_report.pdf"));
        }
    }

    @Test
    public void canGenerateLowDNAYieldReport() throws DRException, IOException {
        JasperReportBuilder report = generateNotAnalysableCPCTReport(0.6, NotAnalysableReason.LOW_DNA_YIELD);
        assertNotNull(report);

        if (WRITE_TO_PDF) {
            report.toPdf(new FileOutputStream(REPORT_BASE_DIR + File.separator + "hmf_low_dna_yield_report.pdf"));
        }
    }

    @Test
    public void canGeneratePostDNAIsolationFailReport() throws DRException, IOException {
        JasperReportBuilder report = generateNotAnalysableCPCTReport(0.6, NotAnalysableReason.POST_ANALYSIS_FAIL);
        assertNotNull(report);

        if (WRITE_TO_PDF) {
            report.toPdf(new FileOutputStream(REPORT_BASE_DIR + File.separator + "hmf_post_dna_isolation_fail_report.pdf"));
        }
    }

    @NotNull
    private static JasperReportBuilder generateNotAnalysableCPCTReport(final double pathologyTumorEstimate,
            @NotNull final NotAnalysableReason reason) {
        NotAnalysedPatientReport patientReport = ImmutableNotAnalysedPatientReport.of(testSampleReport(pathologyTumorEstimate),
                reason,
                NotAnalysableStudy.CPCT,
                Optional.empty(),
                testBaseReportData().signaturePath(),
                testBaseReportData().logoRVAPath());

        return PDFWriter.generateNotAnalysableReport(patientReport);
    }
}