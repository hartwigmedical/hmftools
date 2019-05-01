package com.hartwig.hmftools.patientreporter.cfreport;

import com.hartwig.hmftools.common.ecrf.projections.ImmutablePatientTumorLocation;
import com.hartwig.hmftools.patientreporter.*;

import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailStudy;
import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.Locale;
import java.util.Optional;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testBaseReportData;
import static org.junit.Assert.assertNotNull;

public class CFReportWriterTest {

    private static final boolean WRITE_TO_PDF = true;
    private static final boolean TIMESTAMP_FILES = true;

    private static final String REPORT_BASE_DIR = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";
    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("dd-MMM-yyyy", Locale.ENGLISH);

    @Test
    public void canGeneratePatientReportForCOLO829() throws IOException {
        String filename = WRITE_TO_PDF ? getReportFilePath("hmf_test_sequence_report.pdf") : null;

        AnalysedPatientReport colo829Report = ExampleAnalysisTestFactory.buildCOLO829();
        new CFReportWriter().writeAnalysedPatientReport(colo829Report, filename);
    }

    @Test
    public void canGeneratePatientReportForCompletelyFilledInReport() throws IOException {
        String filename = WRITE_TO_PDF ? getReportFilePath("hmf_full_test_sequence_report.pdf") : null;

        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledIn();
        new CFReportWriter().writeAnalysedPatientReport(patientReport, filename);
    }

    @Test
    public void canGenerateLowTumorPercentageReport() throws IOException {
        String filename = WRITE_TO_PDF ? getReportFilePath("hmf_low_tumor_percentage_report.pdf") : null;
        generateQCFailCPCTReport(0.1, null, QCFailReason.LOW_TUMOR_PERCENTAGE, filename);
    }

    @Test
    public void canGenerateLowDNAYieldReport() throws IOException {
        String filename = WRITE_TO_PDF ? getReportFilePath("hmf_low_dna_yield_report.pdf") : null;
        generateQCFailCPCTReport(0.6, null, QCFailReason.LOW_DNA_YIELD, filename);
    }

    @Test
    public void canGeneratePostDNAIsolationFailReport() throws IOException {
        String filename = WRITE_TO_PDF ? getReportFilePath("hmf_post_dna_isolation_fail_report.pdf") : null;
        generateQCFailCPCTReport(0.6, null, QCFailReason.POST_ANALYSIS_FAIL, filename);

    }

    @Test
    public void canGenerateLowMolecularTumorPercentage() throws IOException {
        String filename = WRITE_TO_PDF ? getReportFilePath("hmf_low_molecular_tumor_percentage_report.pdf") : null;
        generateQCFailCPCTReport(null, 0.15, QCFailReason.SHALLOW_SEQ_LOW_PURITY, filename);
    }

    @NotNull
    private static void generateQCFailCPCTReport(@Nullable Double pathologyTumorPercentage,
                                                                @Nullable Double shallowSeqPurity, @NotNull QCFailReason reason, @Nullable String filename) throws IOException {
        SampleReport sampleReport = ImmutableSampleReport.of("CPCT02991111T",
                "A1",
                "A2",
                ImmutablePatientTumorLocation.of("CPCT02991111", "Skin", "Melanoma"),
                shallowSeqPurity != null ? PatientReportFormat.formatPercent(shallowSeqPurity) : "not determined",
                pathologyTumorPercentage != null ? PatientReportFormat.formatPercent(pathologyTumorPercentage) : "not determined",
                LocalDate.parse("05-Jan-2018", DATE_FORMATTER),
                LocalDate.parse("01-Jan-2018", DATE_FORMATTER),
                "PREP013V23-QC037V20-SEQ008V25",
                "HMF Testing Center",
                "COLO-001-002",
                "ContactMe",
                "contact@me.com",
                "ABC",
                "123456");

        QCFailReport patientReport = ImmutableQCFailReport.of(sampleReport,
                reason,
                QCFailStudy.CPCT,
                Optional.empty(),
                testBaseReportData().signaturePath(),
                testBaseReportData().logoRVAPath());

        CFReportWriter reportWriter = new CFReportWriter();
        reportWriter.writeQCFailReport(patientReport, filename);

    }

    @NotNull
    private static String getReportFilePath(@NotNull String filename) {

        if (TIMESTAMP_FILES) {
            int extensionStart = filename.lastIndexOf('.');
            filename = filename.substring(0, extensionStart) + "_" + System.currentTimeMillis() + filename.substring(extensionStart);
        }

        return REPORT_BASE_DIR + File.separator + filename;

    }

}
