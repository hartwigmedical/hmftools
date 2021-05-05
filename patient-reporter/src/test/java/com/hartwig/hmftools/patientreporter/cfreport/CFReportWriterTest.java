package com.hartwig.hmftools.patientreporter.cfreport;

import java.io.File;
import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.Locale;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.clinical.ImmutablePatientPrimaryTumor;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;
import com.hartwig.hmftools.common.lims.hospital.HospitalContactData;
import com.hartwig.hmftools.common.lims.hospital.ImmutableHospitalContactData;
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
import com.hartwig.hmftools.patientreporter.germline.GermlineReportingModel;
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

    @Test
    public void canGeneratePatientReportForCOLO829() throws IOException {
        AnalysedPatientReport colo829Report = ExampleAnalysisTestFactory.buildCOLO829("PNT00012345T",
                false,
                COLO_COMMENT_STRING,
                PatientReporterTestFactory.createCOLOCohortConfig());

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(colo829Report, testReportFilePath(colo829Report));
    }

    @Test
    public void canGeneratePatientReportForCOLO829Corrected() throws IOException {
        AnalysedPatientReport colo829Report = ExampleAnalysisTestFactory.buildCOLO829("PNT00012345T",
                true,
                COLO_COMMENT_STRING_CORRECTED,
                PatientReporterTestFactory.createCOLOCohortConfig());

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(colo829Report, testReportFilePath(colo829Report));
    }

    @Test
    public void canGeneratePatientReportForCOLO829WithGermline() throws IOException {
        AnalysedPatientReport colo829Report = ExampleAnalysisTestFactory.buildWithCOLO829Data("PNT00012345T_GERMLINE",
                false,
                COLO_COMMENT_STRING,
                QsFormNumber.FOR_080.display(),
                true,
                1D,
                true,
                true,
                PatientReporterTestFactory.createCOLOCohortConfig());

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(colo829Report, testReportFilePath(colo829Report));
    }

    @Test
    public void canGeneratePatientReportForCOLO829BelowDetectionThreshold() throws IOException {
        AnalysedPatientReport colo829Report = ExampleAnalysisTestFactory.buildWithCOLO829Data("PNT00012345T_NO_TUMOR",
                false,
                COLO_COMMENT_STRING,
                QsFormNumber.FOR_209.display(),
                false,
                0.23,
                false,
                false,
                PatientReporterTestFactory.createCOLOCohortConfig());

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(colo829Report, testReportFilePath(colo829Report));
    }

    @Test
    public void canGeneratePatientReportForCOLO829InsufficientTCP() throws IOException {
        AnalysedPatientReport colo829Report = ExampleAnalysisTestFactory.buildWithCOLO829Data("PNT00012345T_INSUFFICIENT_TUMOR",
                false,
                COLO_COMMENT_STRING,
                QsFormNumber.FOR_209.display(),
                true,
                0.19,
                false,
                false,
                PatientReporterTestFactory.createCOLOCohortConfig());

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(colo829Report, testReportFilePath(colo829Report));
    }

    @Test
    public void canGeneratePatientReportForCPCTSample() throws IOException {
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledInAndReliablePurity("CPCT01_FULL",
                FULL_TABLES_COMMENT_STRING,
                PatientReporterTestFactory.createCPCTCohortConfig());

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
    }

    @Test
    public void canGeneratePatientReportForCORESample() throws IOException {
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledInAndReliablePurity("CORE01_FULL",
                null,
                PatientReporterTestFactory.createCORECohortConfig());

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
    }

    @Test
    public void canGeneratePatientReportForWIDESample() throws IOException {
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledInAndReliablePurity("WIDE01_FULL",
                FULL_TABLES_COMMENT_STRING,
                PatientReporterTestFactory.createWIDECohortConfig());

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
    }

    @Test
    public void canGeneratePatientReportForCOREDBSample() throws IOException {
        AnalysedPatientReport patientReport =
                ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledInAndReliablePurity("COREDB01_FULL",
                        FULL_TABLES_COMMENT_STRING,
                        PatientReporterTestFactory.createCOREDBCohortConfig());

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
    }

    @Test
    public void canGeneratePatientReportForBelowDetectionSample() throws IOException {
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledIn("CPCT01_NO_TUMOR_FOR-209",
                FULL_TABLES_COMMENT_STRING,
                false,
                1D,
                PatientReporterTestFactory.createCPCTCohortConfig());

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
    }

    @Test
    public void canGeneratePatientReportForInsufficientTCPSample() throws IOException {
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledIn(
                "CPCT01_INSUFFICIENT_TUMOR-FOR-209",
                FULL_TABLES_COMMENT_STRING,
                true,
                0.19,
                PatientReporterTestFactory.createCPCTCohortConfig());

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
    }

    @Test
    public void canGenerateInsufficientDNAReport() throws IOException {
        generateQCFailCPCTReport("CPCT01_insufficient_dna-FOR-082",
                Lims.NOT_PERFORMED_STRING,
                null,
                QCFailReason.INSUFFICIENT_DNA,
                false,
                COMMENT_STRING_QC_FAIL,
                PatientReporterTestFactory.createCPCTCohortConfig());
    }

    @Test
    public void canGenerateCorrectedInsufficientDNAReport() throws IOException {
        generateQCFailCPCTReport("CPCT01",
                Lims.NOT_PERFORMED_STRING,
                null,
                QCFailReason.INSUFFICIENT_DNA,
                true,
                COMMENT_STRING_QC_FAIL_CORRECTED,
                PatientReporterTestFactory.createCPCTCohortConfig());
    }

    @Test
    public void canGenerateCorrectedInsufficientDNAReportCOREDB() throws IOException {
        generateQCFailCPCTReport("COREDB",
                Lims.NOT_PERFORMED_STRING,
                null,
                QCFailReason.INSUFFICIENT_DNA,
                true,
                COMMENT_STRING_QC_FAIL_CORRECTED,
                PatientReporterTestFactory.createCOREDBCohortConfig());
    }

    @Test
    public void canGenerateTechnicalFailureReport() throws IOException {
        generateQCFailCPCTReport("CPCT02-technical_failure-FOR-102",
                "60%",
                null,
                QCFailReason.TECHNICAL_FAILURE,
                false,
                COMMENT_STRING_QC_FAIL,
                PatientReporterTestFactory.createCPCTCohortConfig());
    }

    @Test
    public void canGenerateSufficientTCPQCFailReport() throws IOException {
        generateQCFailCPCTReport("CPCT03-sufficient_tcp_qc_failure-FOR-083",
                "60%",
                "70%",
                QCFailReason.SUFFICIENT_TCP_QC_FAILURE,
                false,
                COMMENT_STRING_QC_FAIL,
                PatientReporterTestFactory.createCPCTCohortConfig());
    }

    @Test
    public void canGenerateInsufficientTCPAfterDeepWGSReport() throws IOException {
        generateQCFailCPCTReport("CPCT04-insufficient_tcp_deep_wgs-FOR-100",
                "22%",
                "18%",
                QCFailReason.INSUFFICIENT_TCP_DEEP_WGS,
                false,
                COMMENT_STRING_QC_FAIL,
                PatientReporterTestFactory.createCPCTCohortConfig());
    }

    @Test
    public void canGenerateInsufficientTCPAfterShallowReport() throws IOException {
        generateQCFailCPCTReport("CPCT05-insufficient_tcp_shallow_wgs-FOR-100",
                "15%",
                null,
                QCFailReason.INSUFFICIENT_TCP_SHALLOW_WGS,
                false,
                COMMENT_STRING_QC_FAIL,
                PatientReporterTestFactory.createCPCTCohortConfig());
    }

    @Test
    public void canGenerateInsufficientTCPAfterShallowReportCORE() throws IOException {
        generateQCFailCPCTReport("CORE01",
                "15%",
                null,
                QCFailReason.INSUFFICIENT_TCP_SHALLOW_WGS,
                false,
                COMMENT_STRING_QC_FAIL,
                PatientReporterTestFactory.createCORECohortConfig());
    }

    @Test
    public void canGenerateInsufficientTCPAfterShallowReportWIDE() throws IOException {
        generateQCFailCPCTReport("WIDE01",
                "15%",
                null,
                QCFailReason.INSUFFICIENT_TCP_SHALLOW_WGS,
                false,
                COMMENT_STRING_QC_FAIL,
                PatientReporterTestFactory.createWIDECohortConfig());
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

    private static void generateQCFailCPCTReport(@NotNull String sampleId, @NotNull String shallowSeqPurity,
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
                .germlineReportingLevel(LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION)
                .reportViralInsertions(true)
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
                .build();

        String filename = testReportFilePath(patientReport);

        CFReportWriter writer = testCFReportWriter();
        writer.writeQCFailReport(patientReport, filename);
    }

    @NotNull
    private static CFReportWriter testCFReportWriter() {
        GermlineReportingModel emptyGermlineReportingModel = new GermlineReportingModel(Lists.newArrayList());
        return new CFReportWriter(WRITE_TO_PDF, emptyGermlineReportingModel);
    }

    @NotNull
    private static String testReportFilePath(@NotNull PatientReport patientReport) {
        String fileName = OutputFileUtil.generateOutputFileNameForReport(patientReport);
        String newFileName = fileName;
        if (TIMESTAMP_FILES) {
            int extensionStart = fileName.lastIndexOf('.');
            newFileName = fileName.substring(0, extensionStart) + "_" + System.currentTimeMillis() + fileName.substring(extensionStart);
        }
        return REPORT_BASE_DIR + File.separator + newFileName;
    }
}
