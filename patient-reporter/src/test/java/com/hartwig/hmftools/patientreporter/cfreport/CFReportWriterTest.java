package com.hartwig.hmftools.patientreporter.cfreport;

import java.io.File;
import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.Locale;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.clinical.ImmutablePatientPrimaryTumor;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.lims.cohort.ImmutableLimsCohortConfigData;
import com.hartwig.hmftools.common.lims.cohort.ImmutableLimsCohortModel;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfigData;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortModel;
import com.hartwig.hmftools.common.lims.hospital.HospitalContactData;
import com.hartwig.hmftools.common.lims.hospital.ImmutableHospitalContactData;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
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
import com.hartwig.hmftools.patientreporter.qcfail.ImmutableQCFailReport;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReport;
import com.hartwig.hmftools.protect.variants.germline.GermlineReportingModel;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Ignore;
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
        AnalysedPatientReport colo829Report =
                ExampleAnalysisTestFactory.buildCOLO829("PNT00012345T", false, COLO_COMMENT_STRING);

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(colo829Report, testReportFilePath(colo829Report));
    }

    @Test
    public void canGeneratePatientReportForCOLO829Corrected() throws IOException {
        AnalysedPatientReport colo829Report =
                ExampleAnalysisTestFactory.buildCOLO829("PNT00012345T", true, COLO_COMMENT_STRING_CORRECTED);

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
                "CPCT");

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
                "CPCT");

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
                "CPCT");

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(colo829Report, testReportFilePath(colo829Report));
    }

    @Test
    public void canGeneratePatientReportForCPCTSample() throws IOException {
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledInAndReliablePurity("CPCT01_FULL",
                FULL_TABLES_COMMENT_STRING,
                "CPCT");

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
    }

    @Test
    public void canGeneratePatientReportForCORESample() throws IOException {
        AnalysedPatientReport patientReport =
                ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledInAndReliablePurity("CORE01_FULL", null, "CORE");

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
    }

    @Test
    public void canGeneratePatientReportForWIDESample() throws IOException {
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledInAndReliablePurity("WIDE01_FULL",
                FULL_TABLES_COMMENT_STRING,
                "WIDE");

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
    }

    @Test
    public void canGeneratePatientReportForCOREDBSample() throws IOException {
        AnalysedPatientReport patientReport =
                ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledInAndReliablePurity("COREDB01_FULL",
                        FULL_TABLES_COMMENT_STRING,
                        "COREDB");

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
    }

    @Test
    public void canGeneratePatientReportForBelowDetectionSample() throws IOException {
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledIn("CPCT01_NO_TUMOR",
                FULL_TABLES_COMMENT_STRING,
                false,
                1D,
                "CPCT");

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
    }

    @Test
    public void canGeneratePatientReportForInsufficientTCPSample() throws IOException {
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledIn("CPCT01_INSUFFICIENT_TUMOR",
                FULL_TABLES_COMMENT_STRING,
                true,
                0.19,
                "CPCT");

        CFReportWriter writer = testCFReportWriter();
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
    }

    @Test
    public void canGenerateInsufficientDNAReport() throws IOException {
        generateQCFailCPCTReport("CPCT01",
                Lims.NOT_PERFORMED_STRING,
                null,
                QCFailReason.INSUFFICIENT_DNA,
                false,
                COMMENT_STRING_QC_FAIL,
                "CPCT");
    }

    @Test
    public void canGenerateCorrectedInsufficientDNAReport() throws IOException {
        generateQCFailCPCTReport("CPCT01",
                Lims.NOT_PERFORMED_STRING,
                null,
                QCFailReason.INSUFFICIENT_DNA,
                true,
                COMMENT_STRING_QC_FAIL_CORRECTED,
                "CPCT");
    }

    @Test
    public void canGenerateCorrectedInsufficientDNAReportCOREDB() throws IOException {
        generateQCFailCPCTReport("COREDB",
                Lims.NOT_PERFORMED_STRING,
                null,
                QCFailReason.INSUFFICIENT_DNA,
                true,
                COMMENT_STRING_QC_FAIL_CORRECTED,
                "COREDB");
    }

    @Test
    public void canGenerateTechnicalFailureReport() throws IOException {
        generateQCFailCPCTReport("CPCT02", "60%", null, QCFailReason.TECHNICAL_FAILURE, false, COMMENT_STRING_QC_FAIL, "CPCT");
    }

    @Test
    public void canGenerateSufficientTCPQCFailReport() throws IOException {
        generateQCFailCPCTReport("CPCT03",
                "60%",
                "70%",
                QCFailReason.SUFFICIENT_TCP_QC_FAILURE,
                false,
                COMMENT_STRING_QC_FAIL,
                "CPCT");
    }

    @Test
    public void canGenerateInsufficientTCPAfterDeepWGSReport() throws IOException {
        generateQCFailCPCTReport("CPCT04",
                "22%",
                "18%",
                QCFailReason.INSUFFICIENT_TCP_DEEP_WGS,
                false,
                COMMENT_STRING_QC_FAIL,
                "CPCT");
    }

    @Test
    public void canGenerateInsufficientTCPAfterShallowReport() throws IOException {
        generateQCFailCPCTReport("CPCT05",
                "15%",
                null,
                QCFailReason.INSUFFICIENT_TCP_SHALLOW_WGS,
                false,
                COMMENT_STRING_QC_FAIL,
                "CPCT");
    }

    @Test
    public void canGenerateInsufficientTCPAfterShallowReportCORE() throws IOException {
        generateQCFailCPCTReport("CORE01",
                "15%",
                null,
                QCFailReason.INSUFFICIENT_TCP_SHALLOW_WGS,
                false,
                COMMENT_STRING_QC_FAIL,
                "CORE");
    }

    @Test
    public void canGenerateInsufficientTCPAfterShallowReportWIDE() throws IOException {
        generateQCFailCPCTReport("WIDE01",
                "15%",
                null,
                QCFailReason.INSUFFICIENT_TCP_SHALLOW_WGS,
                false,
                COMMENT_STRING_QC_FAIL,
                "WIDE");
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
            @NotNull String cohort) throws IOException {
        SampleMetadata sampleMetadata = ImmutableSampleMetadata.builder()
                .patientId("patient")
                .refSampleId("x")
                .refSampleBarcode("FR12123488")
                .tumorSampleId(sampleId)
                .tumorSampleBarcode("FR12345678")
                .build();
        LimsCohortModel cohortConfig = buildTestCohortModel(cohort);

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
                .cohort(cohortConfig.queryCohortData(cohort, sampleId))
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

    @NotNull
    private static LimsCohortModel buildTestCohortModel(@NotNull String cohortString) {
        Map<String, LimsCohortConfigData> cohortData = Maps.newHashMap();
        LimsCohortConfigData config = ImmutableLimsCohortConfigData.builder()
                .cohortId(cohortString)
                .hospitalId(true)
                .reportGermline(false)
                .reportGermlineFlag(false)
                .reportConclusion(false)
                .reportViral(false)
                .requireHospitalId(false)
                .requireHospitalPAId(false)
                .hospitalPersonsStudy(true)
                .hospitalPersonsRequester(false)
                .outputFile(false)
                .submission(false)
                .sidePanelInfo(false)
                .build();
        cohortData.put(cohortString, config);
        return ImmutableLimsCohortModel.builder().limsCohortMap(cohortData).build();
    }
}
