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
import com.hartwig.hmftools.patientreporter.ExampleAnalysisTestFactory;
import com.hartwig.hmftools.patientreporter.ImmutableQCFailReport;
import com.hartwig.hmftools.patientreporter.ImmutableSampleReport;
import com.hartwig.hmftools.patientreporter.QCFailReport;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailStudy;
import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.exception.DRException;

public class PDFWriterTest {

    private static final boolean WRITE_TO_PDF = false;

    private static final String REPORT_BASE_DIR = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";
    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("dd-MMM-yyyy", Locale.ENGLISH);

    @Test
    public void canGeneratePatientReportForCOLO829() throws DRException, IOException {
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.buildCOLO829();

        JasperReportBuilder report = PDFWriter.generateAnalysedPatientReport(patientReport);
        assertNotNull(report);

        if (WRITE_TO_PDF) {
            report.toPdf(new FileOutputStream(REPORT_BASE_DIR + File.separator + "hmf_test_sequence_report.pdf"));
        }
    }

    @Test
    public void canGeneratePatientReportForCompletelyFilledInReport() throws DRException, IOException {
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledIn();

        JasperReportBuilder report = PDFWriter.generateAnalysedPatientReport(patientReport);
        assertNotNull(report);

        if (WRITE_TO_PDF) {
            report.toPdf(new FileOutputStream(REPORT_BASE_DIR + File.separator + "hmf_full_test_sequence_report.pdf"));
        }
    }

    @Test
    public void canGenerateLowTumorPercentageReport() throws DRException, IOException {
        JasperReportBuilder report = generateQCFailCPCTReport(0.1, null, QCFailReason.LOW_TUMOR_PERCENTAGE);
        assertNotNull(report);

        if (WRITE_TO_PDF) {
            report.toPdf(new FileOutputStream(REPORT_BASE_DIR + File.separator + "hmf_low_tumor_percentage_report.pdf"));
        }
    }

    @Test
    public void canGenerateLowDNAYieldReport() throws DRException, IOException {
        JasperReportBuilder report = generateQCFailCPCTReport(0.6, null, QCFailReason.LOW_DNA_YIELD);
        assertNotNull(report);

        if (WRITE_TO_PDF) {
            report.toPdf(new FileOutputStream(REPORT_BASE_DIR + File.separator + "hmf_low_dna_yield_report.pdf"));
        }
    }

    @Test
    public void canGeneratePostDNAIsolationFailReport() throws DRException, IOException {
        JasperReportBuilder report = generateQCFailCPCTReport(0.6, null, QCFailReason.POST_ANALYSIS_FAIL);
        assertNotNull(report);

        if (WRITE_TO_PDF) {
            report.toPdf(new FileOutputStream(REPORT_BASE_DIR + File.separator + "hmf_post_dna_isolation_fail_report.pdf"));
        }
    }

    @Test
    public void canGenerateLowMolecularTumorPercentage() throws DRException, IOException {
        JasperReportBuilder report = generateQCFailCPCTReport(null, 0.15, QCFailReason.SHALLOW_SEQ_LOW_PURITY);
        assertNotNull(report);

        if (WRITE_TO_PDF) {
            report.toPdf(new FileOutputStream(REPORT_BASE_DIR + File.separator + "hmf_low_molecular_tumor_percentage_report.pdf"));
        }
    }

    @NotNull
    private static JasperReportBuilder generateQCFailCPCTReport(@Nullable Double pathologyTumorPercentage,
            @Nullable Double shallowSeqPurity, @NotNull QCFailReason reason) {
        SampleReport sampleReport = ImmutableSampleReport.builder()
                .sampleId("CPCT02991111T")
                .patientTumorLocation(ImmutablePatientTumorLocation.of("CPCT02991111", "Skin", "Melanoma"))
                .refBarcode("FR12123488")
                .tumorBarcode("FR12345678")
                .tumorArrivalDate(LocalDate.parse("05-Jan-2018", DATE_FORMATTER))
                .purityShallowSeq(shallowSeqPurity != null ? PatientReportFormat.formatPercent(shallowSeqPurity) : "not determined")
                .pathologyTumorPercentage(
                        pathologyTumorPercentage != null ? PatientReportFormat.formatPercent(pathologyTumorPercentage) : "not determined")
                .labProcedures("PREP013V23-QC037V20-SEQ008V25")
                .addressee("HMF Testing Center")
                .projectName("COLO-001-002")
                .requesterName("ContactMe")
                .requesterEmail("contact@me.com")
                .submissionId("ABC")
                .hospitalPatientId("123456")
                .hospitalPathologySampleId("A")
                .build();

        QCFailReport patientReport = ImmutableQCFailReport.of(sampleReport,
                reason,
                QCFailStudy.CPCT,
                Optional.empty(),
                testBaseReportData().signaturePath(),
                testBaseReportData().logoRVAPath());

        return PDFWriter.generateQCFailReport(patientReport);
    }
}