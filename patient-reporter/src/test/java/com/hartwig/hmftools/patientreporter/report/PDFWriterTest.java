package com.hartwig.hmftools.patientreporter.report;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testBaseReportData;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.Locale;
import java.util.Optional;

import com.hartwig.hmftools.common.ecrf.projections.ImmutablePatientTumorLocation;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.ImmutableNotAnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.ImmutableSampleReport;
import com.hartwig.hmftools.patientreporter.NotAnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.qcfail.NotAnalysableReason;
import com.hartwig.hmftools.patientreporter.qcfail.NotAnalysableStudy;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.exception.DRException;

public class PDFWriterTest {

    private static final boolean WRITE_TO_PDF = true;

    private static final String REPORT_BASE_DIR = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";
    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("dd-MMM-yyyy", Locale.ENGLISH);

    @Test
    public void canGenerateSequenceReportForCOLO829() throws DRException, IOException {
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.buildCOLO829();

        JasperReportBuilder report = PDFWriter.generatePatientReport(patientReport);
        assertNotNull(report);

        if (WRITE_TO_PDF) {
            report.toPdf(new FileOutputStream(REPORT_BASE_DIR + File.separator + "hmf_test_sequence_report.pdf"));
        }
    }

    @Test
    public void canGenerateSequenceReportForCompletelyFilledInReport() throws DRException, IOException {
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledIn();

        JasperReportBuilder report = PDFWriter.generatePatientReport(patientReport);
        assertNotNull(report);

        if (WRITE_TO_PDF) {
            report.toPdf(new FileOutputStream(REPORT_BASE_DIR + File.separator + "hmf_test_sequence_report.pdf"));
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
    private static JasperReportBuilder generateNotAnalysableCPCTReport(double pathologyTumorEstimate, @NotNull NotAnalysableReason reason) {
        SampleReport sampleReport = ImmutableSampleReport.of("CPCT02991111T",
                ImmutablePatientTumorLocation.of("CPCT02991111", "Skin", "Melanoma"),
                pathologyTumorEstimate,
                LocalDate.parse("05-Jan-2018", DATE_FORMATTER),
                LocalDate.parse("01-Jan-2018", DATE_FORMATTER),
                "PREP013V23-QC037V20-SEQ008V25",
                "HMF Testing Center");

        NotAnalysedPatientReport patientReport = ImmutableNotAnalysedPatientReport.of(sampleReport,
                reason,
                NotAnalysableStudy.CPCT,
                Optional.empty(),
                testBaseReportData().signaturePath(),
                testBaseReportData().logoRVAPath());

        return PDFWriter.generateNotAnalysableReport(patientReport);
    }
}