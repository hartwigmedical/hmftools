package com.hartwig.hmftools.patientreporter.cfreport;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testReportData;

import java.io.File;
import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.Locale;

import com.hartwig.hmftools.common.ecrf.projections.ImmutablePatientTumorLocation;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.hospital.HospitalContactData;
import com.hartwig.hmftools.common.lims.hospital.ImmutableHospitalContactData;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.ExampleAnalysisTestFactory;
import com.hartwig.hmftools.patientreporter.ImmutableSampleMetadata;
import com.hartwig.hmftools.patientreporter.ImmutableSampleReport;
import com.hartwig.hmftools.patientreporter.SampleMetadata;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.qcfail.ImmutableQCFailReport;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReport;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class CFReportWriterTest {

    private static final boolean WRITE_TO_PDF = false;
    private static final boolean TIMESTAMP_FILES = false;

    private static final String REPORT_BASE_DIR = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";
    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("dd-MMM-yyyy", Locale.ENGLISH);
    private static final String COMMENT_STRING = "This is a test report and is based off COLO829";
    private static final String COMMENT_STRING_CORRECTED = "This is a test corrected report and is based off COLO829";

    private static final String COMMENT_STRING_FAIL = "This is a test qc fail report";
    private static final String COMMENT_STRING_FAIL_CORRECTED = "This is a test corrected qc fail report";

    @Test
    public void canGeneratePatientReportForCOLO829() throws IOException {
        AnalysedPatientReport colo829Report = ExampleAnalysisTestFactory.buildCOLO829(false, COMMENT_STRING);

        CFReportWriter writer = new CFReportWriter(WRITE_TO_PDF);
        writer.writeAnalysedPatientReport(colo829Report, testReportFilePath("hmf_colo829_sequence_report.pdf"));
    }

    @Test
    public void canGeneratePatientReportForCOLO829Corrected() throws IOException {
        AnalysedPatientReport colo829Report = ExampleAnalysisTestFactory.buildCOLO829(true, COMMENT_STRING_CORRECTED);

        CFReportWriter writer = new CFReportWriter(WRITE_TO_PDF);
        writer.writeAnalysedPatientReport(colo829Report, testReportFilePath("hmf_colo829_sequence_report_corrected.pdf"));
    }

    @Test
    public void canGeneratePatientReportForCPCTSample() throws IOException {
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledIn("CPCT01000001T");

        CFReportWriter writer = new CFReportWriter(WRITE_TO_PDF);
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath("hmf_cpct_sequence_report.pdf"));
    }

    @Test
    public void canGeneratePatientReportForCORESample() throws IOException {
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledIn("CORE01000001T");

        CFReportWriter writer = new CFReportWriter(WRITE_TO_PDF);
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath("hmf_core_sequence_report.pdf"));
    }

    @Test
    public void canGeneratePatientReportForWIDESample() throws IOException {
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledIn("WIDE01000001T");

        CFReportWriter writer = new CFReportWriter(WRITE_TO_PDF);
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath("hmf_wide_sequence_report.pdf"));
    }

    @Test
    public void canGeneratePatientReportForBelowDetectionSample() throws IOException {
        AnalysedPatientReport patientReport =
                ExampleAnalysisTestFactory.buildAnalysisWithAllTablesForBelowDetectionLimitSample("CPCT01000001T");

        CFReportWriter writer = new CFReportWriter(WRITE_TO_PDF);
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath("hmf_below_detection_limit_sequence_report.pdf"));
    }

    @Test
    public void canGenerateLowDNAYieldReport() throws IOException {
        generateQCFailCPCTReport("CPCT01000001T",
                "60%",
                QCFailReason.LOW_DNA_YIELD,
                testReportFilePath("hmf_low_dna_yield_report.pdf"),
                false,
                COMMENT_STRING_FAIL);
    }

    @Test
    public void canGenerateLowDNAYieldReportCorrected() throws IOException {
        generateQCFailCPCTReport("CPCT01000001T",
                "60%",
                QCFailReason.LOW_DNA_YIELD,
                testReportFilePath("hmf_low_dna_yield_report_corrected.pdf"),
                true,
                COMMENT_STRING_FAIL_CORRECTED);
    }

    @Test
    public void canGenerateBelowDetectionThresholdWithoutGenomicAlterationsReport() throws IOException {
        generateQCFailCPCTReport("CPCT01000001T",
                "60%",
                QCFailReason.BELOW_DETECTION_THRESHOLD,
                testReportFilePath("hmf_below_detection_without_genomic_alteration_report.pdf"),
                false,
                COMMENT_STRING_FAIL);
    }

    @Test
    public void canGenerateLabFailureReport() throws IOException {
        generateQCFailCPCTReport("CPCT01000001T",
                "60%",
                QCFailReason.LAB_FAILURE,
                testReportFilePath("hmf_lab_failure_report.pdf"),
                false,
                COMMENT_STRING_FAIL);
    }

    @Test
    public void canGenerateInsufficientTissue() throws IOException {
        generateQCFailCPCTReport("CPCT01000001T",
                "60%",
                QCFailReason.INSUFFICIENT_TISSUE,
                testReportFilePath("hmf_insufficient_tissue_report.pdf"),
                false,
                COMMENT_STRING_FAIL);
    }

    @Test
    public void canGeneratePostDNAIsolationFailReport() throws IOException {
        generateQCFailCPCTReport("CPCT01000001T",
                "60%",
                QCFailReason.POST_ANALYSIS_FAIL,
                testReportFilePath("hmf_post_dna_isolation_fail_report.pdf"),
                false,
                COMMENT_STRING_FAIL);
    }

    @Test
    public void canGenerateLowMolecularTumorPercentageCORE() throws IOException {
        generateQCFailCPCTReport("CORE01000001T",
                "15%",
                QCFailReason.SHALLOW_SEQ_LOW_PURITY,
                testReportFilePath("hmf_low_molecular_tumor_percentage_core_report.pdf"),
                false,
                COMMENT_STRING_FAIL);
    }

    @Test
    public void canGenerateLowMolecularTumorPercentageWIDE() throws IOException {
        generateQCFailCPCTReport("WIDE01000001T",
                "15%",
                QCFailReason.SHALLOW_SEQ_LOW_PURITY,
                testReportFilePath("hmf_low_molecular_tumor_percentage_wide_report.pdf"),
                false,
                COMMENT_STRING_FAIL);
    }

    @Test
    public void canGenerateLowMolecularTumorPercentageCPCT() throws IOException {
        generateQCFailCPCTReport("CPCT01000001T",
                "15%",
                QCFailReason.SHALLOW_SEQ_LOW_PURITY,
                testReportFilePath("hmf_low_molecular_tumor_percentage_cpct_report.pdf"),
                false,
                COMMENT_STRING_FAIL);
    }

    @NotNull
    private static HospitalContactData createTestHospitalContactData() {
        return ImmutableHospitalContactData.builder()
                .hospitalPI("AB")
                .requesterName("Paul")
                .requesterEmail("paul@hartwig.com")
                .hospitalName("HMF Testing Center")
                .hospitalAddress("1000 AB AMSTERDAM")
                .build();
    }

    private static void generateQCFailCPCTReport(@NotNull String sampleId, @Nullable String shallowSeqPurity, @NotNull QCFailReason reason,
            @NotNull String filename, boolean correctedReport, @NotNull String comments) throws IOException {
        SampleMetadata sampleMetadata = ImmutableSampleMetadata.builder()
                .refSampleId("x")
                .refSampleBarcode("FR12123488")
                .tumorSampleId(sampleId)
                .tumorSampleBarcode("FR12345678")
                .build();

        SampleReport sampleReport = ImmutableSampleReport.builder()
                .sampleMetadata(sampleMetadata)
                .patientTumorLocation(ImmutablePatientTumorLocation.of(Strings.EMPTY, "Skin", "Melanoma"))
                .refArrivalDate(LocalDate.parse("10-Jan-2019", DATE_FORMATTER))
                .tumorArrivalDate(LocalDate.parse("05-Jan-2019", DATE_FORMATTER))
                .shallowSeqPurityString(shallowSeqPurity != null ? shallowSeqPurity : Lims.NOT_PERFORMED_STRING)
                .labProcedures("PREP013V23-QC037V20-SEQ008V25")
                .cohort("A")
                .projectName("COLO-001-002")
                .submissionId("ABC")
                .hospitalContactData(createTestHospitalContactData())
                .hospitalPatientId("123456")
                .hospitalPathologySampleId("A")
                .build();

        QCFailReport patientReport = ImmutableQCFailReport.builder()
                .sampleReport(sampleReport)
                .reason(reason)
                .wgsPurityString(null)
                .comments(comments)
                .isCorrectedReport(correctedReport)
                .signaturePath(testReportData().signaturePath())
                .logoRVAPath(testReportData().logoRVAPath())
                .logoCompanyPath(testReportData().logoCompanyPath())
                .build();

        CFReportWriter writer = new CFReportWriter(WRITE_TO_PDF);
        writer.writeQCFailReport(patientReport, filename);
    }

    @NotNull
    private static String testReportFilePath(@NotNull String filename) {
        String newFileName = filename;
        if (TIMESTAMP_FILES) {
            int extensionStart = filename.lastIndexOf('.');
            newFileName = filename.substring(0, extensionStart) + "_" + System.currentTimeMillis() + filename.substring(extensionStart);
        }
        return REPORT_BASE_DIR + File.separator + newFileName;
    }
}
