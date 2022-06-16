package com.hartwig.hmftools.patientreporter.cfreport;

import java.io.File;
import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.Locale;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.clinical.ImmutablePatientPrimaryTumor;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortTestFactory;
import com.hartwig.hmftools.common.lims.hospital.HospitalContactData;
import com.hartwig.hmftools.common.lims.hospital.ImmutableHospitalContactData;
import com.hartwig.hmftools.common.peach.ImmutablePeachGenotype;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.utils.DataUtil;
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
import com.hartwig.hmftools.patientreporter.panel.ImmutablePanelFailReport;
import com.hartwig.hmftools.patientreporter.panel.ImmutablePanelReport;
import com.hartwig.hmftools.patientreporter.panel.PanelFailReason;
import com.hartwig.hmftools.patientreporter.panel.PanelFailReport;
import com.hartwig.hmftools.patientreporter.panel.PanelReport;
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
    private static final String COLO_COMMENT_STRING = "This is a test report and is based on COLO829. Where is referred to CKB, "
            + "VICC evidence is listed due to licensing restrictions.";
    private static final String COLO_COMMENT_STRING_CORRECTED = "This is a corrected test report and is based on COLO829";
    private static final String FULL_TABLES_COMMENT_STRING = "This is a test report with all tables filled in";

    private static final String COMMENT_STRING_QC_FAIL = "This is a test QC fail report";
    private static final String COMMENT_STRING_QC_FAIL_CORRECTED = "This is a corrected test QC fail report";

    private static final String UDI_DI = "(01) 8720299486010(8012)v5.25";

    @Test
    public void canGeneratePatientReportForCOLO829() throws IOException {
        ExampleAnalysisConfig config = new ExampleAnalysisConfig.Builder().sampleId("PNT00012345T")
                .comments(COLO_COMMENT_STRING)
                .limsCohortConfig(LimsCohortTestFactory.createCOLOCohortConfig())
                .build();
        AnalysedPatientReport colo829Report = ExampleAnalysisTestFactory.createWithCOLO829Data(config, PurpleQCStatus.FAIL_CONTAMINATION);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(colo829Report, testReportFilePath(colo829Report));
        writer.writeJsonAnalysedFile(colo829Report, REPORT_BASE_DIR);
    }

    @Test
    public void canGeneratePatientReportForCOLO829DisabledConfig() throws IOException {
        ExampleAnalysisConfig config = new ExampleAnalysisConfig.Builder().sampleId("PNT00012345T_disabled_config")
                .comments(COLO_COMMENT_STRING)
                .limsCohortConfig(LimsCohortTestFactory.createAllDisabledCohortConfig("COLO"))
                .build();
        AnalysedPatientReport colo829Report = ExampleAnalysisTestFactory.createWithCOLO829Data(config, PurpleQCStatus.PASS);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(colo829Report, testReportFilePath(colo829Report));
        writer.writeJsonAnalysedFile(colo829Report, REPORT_BASE_DIR);
    }

    @Test
    public void canGeneratePatientReportForCOLO829Corrected() throws IOException {
        ExampleAnalysisConfig config = new ExampleAnalysisConfig.Builder().sampleId("PNT00012345T")
                .isCorrectionReport(true)
                .comments(COLO_COMMENT_STRING_CORRECTED)
                .limsCohortConfig(LimsCohortTestFactory.createCOLOCohortConfig())
                .build();
        AnalysedPatientReport colo829Report = ExampleAnalysisTestFactory.createWithCOLO829Data(config, PurpleQCStatus.PASS);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(colo829Report, testReportFilePath(colo829Report));
        writer.writeJsonAnalysedFile(colo829Report, REPORT_BASE_DIR);
    }

    @Test
    public void canGeneratePatientReportForCOLO829WithGermline() throws IOException {
        ExampleAnalysisConfig config = new ExampleAnalysisConfig.Builder().sampleId("PNT00012345T_GERMLINE")
                .comments(COLO_COMMENT_STRING)
                .reportGermline(true)
                .limsCohortConfig(LimsCohortTestFactory.createCOLOCohortConfig())
                .build();
        AnalysedPatientReport colo829Report = ExampleAnalysisTestFactory.createWithCOLO829Data(config, PurpleQCStatus.PASS);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(colo829Report, testReportFilePath(colo829Report));
        writer.writeJsonAnalysedFile(colo829Report, REPORT_BASE_DIR);
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
        AnalysedPatientReport colo829Report = ExampleAnalysisTestFactory.createWithCOLO829Data(config, PurpleQCStatus.PASS);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(colo829Report, testReportFilePath(colo829Report));
        writer.writeJsonAnalysedFile(colo829Report, REPORT_BASE_DIR);
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
        AnalysedPatientReport colo829Report = ExampleAnalysisTestFactory.createWithCOLO829Data(config, PurpleQCStatus.PASS);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(colo829Report, testReportFilePath(colo829Report));
        writer.writeJsonAnalysedFile(colo829Report, REPORT_BASE_DIR);
    }

    @Test
    public void canGeneratePatientReportForCPCTSample() throws IOException {
        ExampleAnalysisConfig config = new ExampleAnalysisConfig.Builder().sampleId("CPCT01_FULL")
                .comments(FULL_TABLES_COMMENT_STRING)
                .limsCohortConfig(LimsCohortTestFactory.createCPCTCohortConfig())
                .build();
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.createAnalysisWithAllTablesFilledIn(config, PurpleQCStatus.PASS);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
        writer.writeJsonAnalysedFile(patientReport, REPORT_BASE_DIR);
    }

    @Test
    public void canGeneratePatientReportForACTINSample() throws IOException {
        ExampleAnalysisConfig config = new ExampleAnalysisConfig.Builder().sampleId("ACTN_FULL")
                .limsCohortConfig(LimsCohortTestFactory.createACTINCohortConfig())
                .build();
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.createAnalysisWithAllTablesFilledIn(config, PurpleQCStatus.PASS);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
        writer.writeJsonAnalysedFile(patientReport, REPORT_BASE_DIR);
    }

    @Test
    public void canGeneratePatientReportForCORESample() throws IOException {
        ExampleAnalysisConfig config = new ExampleAnalysisConfig.Builder().sampleId("CORE01_FULL")
                .limsCohortConfig(LimsCohortTestFactory.createCORECohortConfig())
                .build();
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.createAnalysisWithAllTablesFilledIn(config, PurpleQCStatus.PASS);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
        writer.writeJsonAnalysedFile(patientReport, REPORT_BASE_DIR);
    }

    @Test
    public void canGeneratePatientReportForWIDESample() throws IOException {
        ExampleAnalysisConfig config = new ExampleAnalysisConfig.Builder().sampleId("WIDE01_FULL")
                .comments(FULL_TABLES_COMMENT_STRING)
                .limsCohortConfig(LimsCohortTestFactory.createWIDECohortConfig())
                .build();
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.createAnalysisWithAllTablesFilledIn(config, PurpleQCStatus.PASS);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
        writer.writeJsonAnalysedFile(patientReport, REPORT_BASE_DIR);
    }

    @Test
    public void canGeneratePatientReportForCOREDBSample() throws IOException {
        ExampleAnalysisConfig config = new ExampleAnalysisConfig.Builder().sampleId("COREDB01_FULL")
                .comments(FULL_TABLES_COMMENT_STRING)
                .limsCohortConfig(LimsCohortTestFactory.createCOREDBCohortConfig())
                .build();
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.createAnalysisWithAllTablesFilledIn(config, PurpleQCStatus.PASS);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
        writer.writeJsonAnalysedFile(patientReport, REPORT_BASE_DIR);
    }

    @Test
    public void canGeneratePatientReportForBelowDetectionSample() throws IOException {
        ExampleAnalysisConfig config = new ExampleAnalysisConfig.Builder().sampleId("CPCT01_NO_TUMOR_FOR-209")
                .comments(FULL_TABLES_COMMENT_STRING)
                .hasReliablePurity(false)
                .qcForNumber(QsFormNumber.FOR_209)
                .limsCohortConfig(LimsCohortTestFactory.createCPCTCohortConfig())
                .build();
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.createAnalysisWithAllTablesFilledIn(config, PurpleQCStatus.PASS);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
        writer.writeJsonAnalysedFile(patientReport, REPORT_BASE_DIR);
    }

    @Test
    public void canGeneratePatientReportForInsufficientTCPSample() throws IOException {
        ExampleAnalysisConfig config = new ExampleAnalysisConfig.Builder().sampleId("CPCT01_INSUFFICIENT_TUMOR-FOR-209")
                .comments(FULL_TABLES_COMMENT_STRING)
                .impliedTumorPurity(0.19)
                .qcForNumber(QsFormNumber.FOR_209)
                .limsCohortConfig(LimsCohortTestFactory.createCPCTCohortConfig())
                .build();
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.createAnalysisWithAllTablesFilledIn(config, PurpleQCStatus.PASS);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
        writer.writeJsonAnalysedFile(patientReport, REPORT_BASE_DIR);
    }

    @Test
    public void canGenerateInsufficientDNAReport() throws IOException {
        generateQCFailReport("CPCT01_insufficient_dna-FOR-082",
                Lims.NOT_PERFORMED_STRING,
                null,
                QCFailReason.INSUFFICIENT_DNA,
                false,
                false,
                COMMENT_STRING_QC_FAIL,
                LimsCohortTestFactory.createCPCTCohortConfig(),
                PurpleQCStatus.PASS);
    }

    @Test
    public void canGenerateCorrectedInsufficientDNAReport() throws IOException {
        generateQCFailReport("CPCT01",
                Lims.NOT_PERFORMED_STRING,
                null,
                QCFailReason.INSUFFICIENT_DNA,
                true,
                true,
                COMMENT_STRING_QC_FAIL_CORRECTED,
                LimsCohortTestFactory.createCPCTCohortConfig(),
                PurpleQCStatus.PASS);
    }

    @Test
    public void canGenerateCorrectedInsufficientDNAReportCOREDB() throws IOException {
        generateQCFailReport("COREDB",
                Lims.NOT_PERFORMED_STRING,
                null,
                QCFailReason.INSUFFICIENT_DNA,
                true,
                false,
                COMMENT_STRING_QC_FAIL_CORRECTED,
                LimsCohortTestFactory.createCOREDBCohortConfig(),
                PurpleQCStatus.PASS);
    }

    @Test
    public void canGenerateTechnicalFailureReport() throws IOException {
        generateQCFailReport("CPCT02-technical_failure-FOR-102",
                "60%",
                null,
                QCFailReason.TECHNICAL_FAILURE,
                false,
                false,
                COMMENT_STRING_QC_FAIL,
                LimsCohortTestFactory.createCPCTCohortConfig(),
                PurpleQCStatus.PASS);
    }

    @Test
    public void canGenerateSufficientTCPQCFailReport() throws IOException {
        generateQCFailReport("CPCT03-sufficient_tcp_qc_failure-FOR-083_Fail",
                "60%",
                "70%",
                QCFailReason.SUFFICIENT_TCP_QC_FAILURE,
                false,
                false,
                COMMENT_STRING_QC_FAIL,
                LimsCohortTestFactory.createCPCTCohortConfig(),
                PurpleQCStatus.FAIL_CONTAMINATION);
    }

    @Test
    public void canGenerateSufficientTCPQCFailReportPASS() throws IOException {
        generateQCFailReport("CPCT03-sufficient_tcp_qc_failure-FOR-083_Pass",
                "60%",
                "70%",
                QCFailReason.SUFFICIENT_TCP_QC_FAILURE,
                false,
                false,
                COMMENT_STRING_QC_FAIL,
                LimsCohortTestFactory.createCPCTCohortConfig(),
                PurpleQCStatus.PASS);
    }

    @Test
    public void canGenerateInsufficientTCPAfterDeepWGSReport() throws IOException {
        generateQCFailReport("CPCT04-insufficient_tcp_deep_wgs-FOR-100",
                "22%",
                "18%",
                QCFailReason.INSUFFICIENT_TCP_DEEP_WGS,
                false,
                false,
                COMMENT_STRING_QC_FAIL,
                LimsCohortTestFactory.createCPCTCohortConfig(),
                PurpleQCStatus.PASS);
    }

    @Test
    public void canGenerateInsufficientTCPAfterShallowReport() throws IOException {
        generateQCFailReport("CPCT05-insufficient_tcp_shallow_wgs-FOR-100",
                "15%",
                null,
                QCFailReason.INSUFFICIENT_TCP_SHALLOW_WGS,
                false,
                false,
                COMMENT_STRING_QC_FAIL,
                LimsCohortTestFactory.createCPCTCohortConfig(),
                PurpleQCStatus.PASS);
    }

    @Test
    public void canGenerateInsufficientTCPAfterShallowReportCORE() throws IOException {
        generateQCFailReport("CORE01",
                "15%",
                null,
                QCFailReason.INSUFFICIENT_TCP_SHALLOW_WGS,
                false,
                false,
                COMMENT_STRING_QC_FAIL,
                LimsCohortTestFactory.createCORECohortConfig(),
                PurpleQCStatus.PASS);
    }

    @Test
    public void canGenerateInsufficientTCPAfterShallowReportWIDE() throws IOException {
        generateQCFailReport("WIDE01",
                "15%",
                null,
                QCFailReason.INSUFFICIENT_TCP_SHALLOW_WGS,
                false,
                false,
                COMMENT_STRING_QC_FAIL,
                LimsCohortTestFactory.createWIDECohortConfig(),
                PurpleQCStatus.PASS);
    }

    @Test
    public void generatePanelReport() throws IOException {
        SampleMetadata sampleMetadata = generateSampleMetadata("Sample_panel");
        ReportData testReportData = PatientReporterTestFactory.loadTestReportData();
        SampleReport sampleReport = generateSampleReport(sampleMetadata);

        PanelReport patientReport = ImmutablePanelReport.builder()
                .sampleReport(sampleReport)
                .qsFormNumber("form")
                .VCFFilename("test.vcf")
                .isCorrectedReport(false)
                .isCorrectedReportExtern(false)
                .signaturePath(testReportData.signaturePath())
                .logoCompanyPath(testReportData.logoCompanyPath())
                .reportDate(DataUtil.formatDate(LocalDate.now()))
                .isWGSreport(false)
                .comments("This is a test report")
                .build();

        String filename = testReportFilePathPanel(patientReport);

        CFReportWriter writer = testCFReportWriter();
        writer.writePanelAnalysedReport(patientReport, filename);
    }

    @Test
    public void generateFailPanelReport() throws IOException {
        SampleMetadata sampleMetadata = generateSampleMetadata("sample_panel_failed");
        ReportData testReportData = PatientReporterTestFactory.loadTestReportData();

        SampleReport sampleReport = generateSampleReport(sampleMetadata);
        PanelFailReport patientReport = ImmutablePanelFailReport.builder()
                .sampleReport(sampleReport)
                .qsFormNumber("form")
                .panelFailReason(PanelFailReason.PANEL_FAILURE)
                .isCorrectedReport(false)
                .isCorrectedReportExtern(false)
                .signaturePath(testReportData.signaturePath())
                .logoCompanyPath(testReportData.logoCompanyPath())
                .reportDate(DataUtil.formatDate(LocalDate.now()))
                .isWGSreport(false)
                .comments("This is a test report")
                .build();

        String filename = testReportFilePathPanel(patientReport);

        CFReportWriter writer = testCFReportWriter();
        writer.writePanelQCFailReport(patientReport, filename);
    }

    private static SampleMetadata generateSampleMetadata(@NotNull String sampleId) {
        return ImmutableSampleMetadata.builder()
                .refSampleId("x")
                .refSampleBarcode("FR12123488")
                .tumorSampleId(sampleId)
                .tumorSampleBarcode("FR12345678")
                .sampleNameForReport(sampleId)
                .build();
    }

    private static SampleReport generateSampleReport(@NotNull SampleMetadata sampleMetadata) {
        return ImmutableSampleReport.builder()
                .sampleMetadata(sampleMetadata)
                .patientPrimaryTumor(ImmutablePatientPrimaryTumor.builder()
                        .patientIdentifier("test")
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
                .shallowSeqPurityString("")
                .labProcedures("PREP013V23-QC037V20-SEQ008V25")
                .cohort(LimsCohortTestFactory.createCOREDBCohortConfig())
                .projectName("TEST-001-002")
                .submissionId("SUBM")
                .hospitalContactData(createTestHospitalContactData())
                .hospitalPatientId("HOSP1")
                .hospitalPathologySampleId("PA1")
                .build();
    }
    private static void generateQCFailReport(@NotNull String sampleId, @NotNull String shallowSeqPurity, @Nullable String wgsPurityString,
            @NotNull QCFailReason reason, boolean correctedReport, boolean correctionReportExtern, @NotNull String comments,
            @NotNull LimsCohortConfig limsCohortConfig, @NotNull PurpleQCStatus purpleQCStatus) throws IOException {
        SampleMetadata sampleMetadata = generateSampleMetadata(sampleId);

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
                .isCorrectedReportExtern(correctionReportExtern)
                .signaturePath(testReportData.signaturePath())
                .logoRVAPath(testReportData.logoRVAPath())
                .logoCompanyPath(testReportData.logoCompanyPath())
                .udiDi(UDI_DI)
                .peachGenotypes(createTestPeachGenotypes())
                .purpleQC(Sets.newHashSet(purpleQCStatus))
                .reportDate(DataUtil.formatDate(LocalDate.now()))
                .isWGSreport(true)
                .build();

        String filename = testReportFilePath(patientReport);

        CFReportWriter writer = testCFReportWriter();
        writer.writeQCFailReport(patientReport, filename);
        writer.writeJsonFailedFile(patientReport, REPORT_BASE_DIR);
    }

    @NotNull
    private static List<PeachGenotype> createTestPeachGenotypes() {
        return Lists.newArrayList(ImmutablePeachGenotype.builder()
                .gene("DPYD")
                .haplotype("*1_HOM")
                .function("Normal Function")
                .linkedDrugs("5-Fluorouracil;Capecitabine;Tegafur")
                .urlPrescriptionInfo("https://www.pharmgkb.org/chemical/PA128406956/guidelineAnnotation/PA166104939;"
                        + "https://www.pharmgkb.org/chemical/PA448771/guidelineAnnotation/PA166104963;"
                        + "https://www.pharmgkb.org/chemical/PA452620/guidelineAnnotation/PA166104944")
                .panelVersion("PGx_min_DPYD_v0.3")
                .repoVersion("1.0")
                .build());
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

    @NotNull
    private static String testReportFilePathPanel(@NotNull com.hartwig.hmftools.patientreporter.PanelReport patientReport) {
        String fileName = OutputFileUtil.generateOutputFileNameForPdfPanelResultReport(patientReport);
        String newFileName = fileName;
        if (TIMESTAMP_FILES) {
            int extensionStart = fileName.lastIndexOf('.');
            newFileName = fileName.substring(0, extensionStart) + "_" + System.currentTimeMillis() + fileName.substring(extensionStart);
        }
        return REPORT_BASE_DIR + File.separator + newFileName;
    }
}
