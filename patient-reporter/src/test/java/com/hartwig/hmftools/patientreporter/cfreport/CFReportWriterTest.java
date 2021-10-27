package com.hartwig.hmftools.patientreporter.cfreport;

import java.io.File;
import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.Locale;

import com.hartwig.hmftools.common.clinical.ImmutablePatientPrimaryTumor;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortTestFactory;
import com.hartwig.hmftools.common.lims.hospital.HospitalContactData;
import com.hartwig.hmftools.common.lims.hospital.ImmutableHospitalContactData;
import com.hartwig.hmftools.patientreporter.ExampleAnalysisConfig;
import com.hartwig.hmftools.patientreporter.ExampleAnalysisTestFactory;
import com.hartwig.hmftools.patientreporter.ImmutableSampleMetadata;
import com.hartwig.hmftools.patientreporter.ImmutableSampleReport;
import com.hartwig.hmftools.patientreporter.OutputFileUtil;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.PatientReporterTestFactory;
import com.hartwig.hmftools.patientreporter.QsFormNumber;
import com.hartwig.hmftools.patientreporter.ReportData;
import com.hartwig.hmftools.patientreporter.SampleMetadata;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.algo.AnalysedPatientReport;
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
    private static final String COLO_COMMENT_STRING = "This is a test report and is based on COLO829";
    private static final String COLO_COMMENT_STRING_CORRECTED = "This is a corrected test report and is based on COLO829";
    private static final String FULL_TABLES_COMMENT_STRING = "This is a test report with all tables filled in";

    private static final String COMMENT_STRING_QC_FAIL = "This is a test QC fail report";
    private static final String COMMENT_STRING_QC_FAIL_CORRECTED = "This is a corrected test QC fail report";

    private static final String UDI_DI = "5.22";

    @Test
    public void canGeneratePatientReportForCOLO829() throws IOException {
        ExampleAnalysisConfig config = new ExampleAnalysisConfig.Builder().sampleId("PNT00012345T")
                .comments(COLO_COMMENT_STRING)
                .limsCohortConfig(LimsCohortTestFactory.createCOLOCohortConfig())
                .build();
        AnalysedPatientReport colo829Report = ExampleAnalysisTestFactory.createWithCOLO829Data(config);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(colo829Report, testReportFilePath(colo829Report));
    }

    @Test
    public void canGeneratePatientReportForCOLO829_disabled_config() throws IOException {
        ExampleAnalysisConfig config = new ExampleAnalysisConfig.Builder().sampleId("PNT00012345T_disables_config")
                .comments(COLO_COMMENT_STRING)
                .limsCohortConfig(LimsCohortTestFactory.createAllDisabledCohortConfig("COLO"))
                .build();
        AnalysedPatientReport colo829Report = ExampleAnalysisTestFactory.createWithCOLO829Data(config);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(colo829Report, testReportFilePath(colo829Report));
    }

    @Test
    public void canGeneratePatientReportForCOLO829Corrected() throws IOException {
        ExampleAnalysisConfig config = new ExampleAnalysisConfig.Builder().sampleId("PNT00012345T")
                .isCorrectionReport(true)
                .comments(COLO_COMMENT_STRING_CORRECTED)
                .limsCohortConfig(LimsCohortTestFactory.createCOLOCohortConfig())
                .build();
        AnalysedPatientReport colo829Report = ExampleAnalysisTestFactory.createWithCOLO829Data(config);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(colo829Report, testReportFilePath(colo829Report));
    }

    @Test
    public void canGeneratePatientReportForCOLO829WithGermline() throws IOException {
        ExampleAnalysisConfig config = new ExampleAnalysisConfig.Builder().sampleId("PNT00012345T_GERMLINE")
                .comments(COLO_COMMENT_STRING)
                .reportGermline(true)
                .limsCohortConfig(LimsCohortTestFactory.createCOLOCohortConfig())
                .build();
        AnalysedPatientReport colo829Report = ExampleAnalysisTestFactory.createWithCOLO829Data(config);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(colo829Report, testReportFilePath(colo829Report));
    }

    @Test
    public void canGeneratePatientReportForCOLO829BelowDetectionThreshold() throws IOException {
        ExampleAnalysisConfig config = new ExampleAnalysisConfig.Builder().sampleId("PNT00012345T_NO_TUMOR")
                .comments(COLO_COMMENT_STRING)
                .qcForNumber(QsFormNumber.FOR_209)
                .hasReliablePurity(false)
                .includeSummary(false)
                .limsCohortConfig(LimsCohortTestFactory.createCOLOCohortConfig())
                .build();
        AnalysedPatientReport colo829Report = ExampleAnalysisTestFactory.createWithCOLO829Data(config);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(colo829Report, testReportFilePath(colo829Report));
    }

    @Test
    public void canGeneratePatientReportForCOLO829InsufficientTCP() throws IOException {
        ExampleAnalysisConfig config = new ExampleAnalysisConfig.Builder().sampleId("PNT00012345T_INSUFFICIENT_TUMOR")
                .comments(COLO_COMMENT_STRING)
                .qcForNumber(QsFormNumber.FOR_209)
                .impliedTumorPurity(0.19)
                .includeSummary(false)
                .hasReliablePurity(true)
                .limsCohortConfig(LimsCohortTestFactory.createCOLOCohortConfig())
                .build();
        AnalysedPatientReport colo829Report = ExampleAnalysisTestFactory.createWithCOLO829Data(config);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(colo829Report, testReportFilePath(colo829Report));
    }

    @Test
    public void canGeneratePatientReportForCPCTSample() throws IOException {
        ExampleAnalysisConfig config = new ExampleAnalysisConfig.Builder().sampleId("CPCT01_FULL")
                .comments(FULL_TABLES_COMMENT_STRING)
                .limsCohortConfig(LimsCohortTestFactory.createCPCTCohortConfig())
                .build();
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.createAnalysisWithAllTablesFilledIn(config);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
    }

    @Test
    public void canGeneratePatientReportForACTINSample() throws IOException {
        ExampleAnalysisConfig config = new ExampleAnalysisConfig.Builder().sampleId("ACTN_FULL")
                .limsCohortConfig(LimsCohortTestFactory.createACTINCohortConfig())
                .build();
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.createAnalysisWithAllTablesFilledIn(config);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
    }

    @Test
    public void canGeneratePatientReportForCORESample() throws IOException {
        ExampleAnalysisConfig config = new ExampleAnalysisConfig.Builder().sampleId("CORE01_FULL")
                .limsCohortConfig(LimsCohortTestFactory.createCORECohortConfig())
                .build();
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.createAnalysisWithAllTablesFilledIn(config);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
    }

    @Test
    public void canGeneratePatientReportForWIDESample() throws IOException {
        ExampleAnalysisConfig config = new ExampleAnalysisConfig.Builder().sampleId("WIDE01_FULL")
                .comments(FULL_TABLES_COMMENT_STRING)
                .limsCohortConfig(LimsCohortTestFactory.createWIDECohortConfig())
                .build();
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.createAnalysisWithAllTablesFilledIn(config);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
    }

    @Test
    public void canGeneratePatientReportForCOREDBSample() throws IOException {
        ExampleAnalysisConfig config = new ExampleAnalysisConfig.Builder().sampleId("COREDB01_FULL")
                .comments(FULL_TABLES_COMMENT_STRING)
                .limsCohortConfig(LimsCohortTestFactory.createCOREDBCohortConfig())
                .build();
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.createAnalysisWithAllTablesFilledIn(config);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
    }

    @Test
    public void canGeneratePatientReportForBelowDetectionSample() throws IOException {
        ExampleAnalysisConfig config = new ExampleAnalysisConfig.Builder().sampleId("CPCT01_NO_TUMOR_FOR-209")
                .comments(FULL_TABLES_COMMENT_STRING)
                .hasReliablePurity(false)
                .qcForNumber(QsFormNumber.FOR_209)
                .limsCohortConfig(LimsCohortTestFactory.createCPCTCohortConfig())
                .build();
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.createAnalysisWithAllTablesFilledIn(config);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
    }

    @Test
    public void canGeneratePatientReportForInsufficientTCPSample() throws IOException {
        ExampleAnalysisConfig config = new ExampleAnalysisConfig.Builder().sampleId("CPCT01_INSUFFICIENT_TUMOR-FOR-209")
                .comments(FULL_TABLES_COMMENT_STRING)
                .impliedTumorPurity(0.19)
                .qcForNumber(QsFormNumber.FOR_209)
                .limsCohortConfig(LimsCohortTestFactory.createCPCTCohortConfig())
                .build();
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.createAnalysisWithAllTablesFilledIn(config);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
    }

    @Test
    public void canGenerateInsufficientDNAReport() throws IOException {
        generateQCFailReport("CPCT01_insufficient_dna-FOR-082",
                Lims.NOT_PERFORMED_STRING,
                null,
                QCFailReason.INSUFFICIENT_DNA,
                false,
                COMMENT_STRING_QC_FAIL,
                LimsCohortTestFactory.createCPCTCohortConfig());
    }

    @Test
    public void canGenerateCorrectedInsufficientDNAReport() throws IOException {
        generateQCFailReport("CPCT01",
                Lims.NOT_PERFORMED_STRING,
                null,
                QCFailReason.INSUFFICIENT_DNA,
                true,
                COMMENT_STRING_QC_FAIL_CORRECTED,
                LimsCohortTestFactory.createCPCTCohortConfig());
    }

    @Test
    public void canGenerateCorrectedInsufficientDNAReportCOREDB() throws IOException {
        generateQCFailReport("COREDB",
                Lims.NOT_PERFORMED_STRING,
                null,
                QCFailReason.INSUFFICIENT_DNA,
                true,
                COMMENT_STRING_QC_FAIL_CORRECTED,
                LimsCohortTestFactory.createCOREDBCohortConfig());
    }

    @Test
    public void canGenerateTechnicalFailureReport() throws IOException {
        generateQCFailReport("CPCT02-technical_failure-FOR-102",
                "60%",
                null,
                QCFailReason.TECHNICAL_FAILURE,
                false,
                COMMENT_STRING_QC_FAIL,
                LimsCohortTestFactory.createCPCTCohortConfig());
    }

    @Test
    public void canGenerateSufficientTCPQCFailReport() throws IOException {
        generateQCFailReport("CPCT03-sufficient_tcp_qc_failure-FOR-083",
                "60%",
                "70%",
                QCFailReason.SUFFICIENT_TCP_QC_FAILURE,
                false,
                COMMENT_STRING_QC_FAIL,
                LimsCohortTestFactory.createCPCTCohortConfig());
    }

    @Test
    public void canGenerateInsufficientTCPAfterDeepWGSReport() throws IOException {
        generateQCFailReport("CPCT04-insufficient_tcp_deep_wgs-FOR-100",
                "22%",
                "18%",
                QCFailReason.INSUFFICIENT_TCP_DEEP_WGS,
                false,
                COMMENT_STRING_QC_FAIL,
                LimsCohortTestFactory.createCPCTCohortConfig());
    }

    @Test
    public void canGenerateInsufficientTCPAfterShallowReport() throws IOException {
        generateQCFailReport("CPCT05-insufficient_tcp_shallow_wgs-FOR-100",
                "15%",
                null,
                QCFailReason.INSUFFICIENT_TCP_SHALLOW_WGS,
                false,
                COMMENT_STRING_QC_FAIL,
                LimsCohortTestFactory.createCPCTCohortConfig());
    }

    @Test
    public void canGenerateInsufficientTCPAfterShallowReportCORE() throws IOException {
        generateQCFailReport("CORE01",
                "15%",
                null,
                QCFailReason.INSUFFICIENT_TCP_SHALLOW_WGS,
                false,
                COMMENT_STRING_QC_FAIL,
                LimsCohortTestFactory.createCORECohortConfig());
    }

    @Test
    public void canGenerateInsufficientTCPAfterShallowReportWIDE() throws IOException {
        generateQCFailReport("WIDE01",
                "15%",
                null,
                QCFailReason.INSUFFICIENT_TCP_SHALLOW_WGS,
                false,
                COMMENT_STRING_QC_FAIL,
                LimsCohortTestFactory.createWIDECohortConfig());
    }

    private static void generateQCFailReport(@NotNull String sampleId, @NotNull String shallowSeqPurity,
            @Nullable String wgsPurityString, @NotNull QCFailReason reason, boolean correctedReport, @NotNull String comments,
            @NotNull LimsCohortConfig limsCohortConfig) throws IOException {
        SampleMetadata sampleMetadata = ImmutableSampleMetadata.builder()
                .refSampleId("x")
                .refSampleBarcode("FR12123488")
                .tumorSampleId(sampleId)
                .tumorSampleBarcode("FR12345678")
                .build();

        SampleReport sampleReport = ImmutableSampleReport.builder()
                .sampleMetadata(sampleMetadata)
                .patientPrimaryTumor(ImmutablePatientPrimaryTumor.builder()
                        .patientIdentifier(sampleId)
                        .location("Skin")
                        .subLocation(Strings.EMPTY)
                        .type("Melanoma")
                        .subType(Strings.EMPTY)
                        .extraDetails(Strings.EMPTY)
                        .isOverridden(false)
                        .build())
                .biopsyLocation(Strings.EMPTY)
                .germlineReportingLevel(LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION)
                .reportViralPresence(true)
                .reportPharmogenetics(true)
                .refArrivalDate(LocalDate.parse("10-Jan-2020", DATE_FORMATTER))
                .tumorArrivalDate(LocalDate.parse("05-Jan-2020", DATE_FORMATTER))
                .shallowSeqPurityString(shallowSeqPurity)
                .labProcedures("PREP013V23-QC037V20-SEQ008V25")
                .cohort(limsCohortConfig)
                .projectName("TEST-001-002")
                .submissionId("SUBM")
                .hospitalContactData(createTestHospitalContactData())
                .hospitalPatientId("HOSP1")
                .hospitalPathologySampleId("PA1")
                .build();

        ReportData testReportData = PatientReporterTestFactory.loadTestReportData();
        QCFailReport patientReport = ImmutableQCFailReport.builder()
                .sampleReport(sampleReport)
                .qsFormNumber(reason.qcFormNumber())
                .reason(reason)
                .wgsPurityString(wgsPurityString)
                .comments(comments)
                .isCorrectedReport(correctedReport)
                .signaturePath(testReportData.signaturePath())
                .logoRVAPath(testReportData.logoRVAPath())
                .logoCompanyPath(testReportData.logoCompanyPath())
                .udiDi(UDI_DI)
                .build();

        String filename = testReportFilePath(patientReport);

        CFReportWriter writer = testCFReportWriter();
        writer.writeQCFailReport(patientReport, filename);
    }

    @NotNull
    private static HospitalContactData createTestHospitalContactData() {
        return ImmutableHospitalContactData.builder()
                .hospitalPI("PI")
                .requesterName("Paul")
                .requesterEmail("paul@hartwig.com")
                .hospitalName("HMF Testing Center")
                .hospitalAddress("1000 AB AMSTERDAM")
                .build();
    }

    @NotNull
    private static CFReportWriter testCFReportWriter() {
        return new CFReportWriter(WRITE_TO_PDF);
    }

    @NotNull
    private static String testReportFilePath(@NotNull PatientReport patientReport) {
        String fileName = OutputFileUtil.generateOutputFileNameForPdfReport(patientReport);
        String newFileName = fileName;
        if (TIMESTAMP_FILES) {
            int extensionStart = fileName.lastIndexOf('.');
            newFileName = fileName.substring(0, extensionStart) + "_" + System.currentTimeMillis() + fileName.substring(extensionStart);
        }
        return REPORT_BASE_DIR + File.separator + newFileName;
    }
}
