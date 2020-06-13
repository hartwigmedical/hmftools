package com.hartwig.hmftools.patientreporter.cfreport;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testReportData;

import java.io.File;
import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.Locale;

import com.hartwig.hmftools.common.ecrf.projections.ImmutablePatientTumorLocation;
import com.hartwig.hmftools.common.lims.hospital.HospitalContactData;
import com.hartwig.hmftools.common.lims.hospital.ImmutableHospitalContactData;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.ExampleAnalysisTestFactory;
import com.hartwig.hmftools.patientreporter.ImmutableSampleMetadata;
import com.hartwig.hmftools.patientreporter.ImmutableSampleReport;
import com.hartwig.hmftools.patientreporter.OutputFileUtil;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.SampleMetadata;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.qcfail.ImmutableQCFailReport;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReport;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CFReportWriterTest {

    private static final boolean WRITE_TO_PDF = false;
    private static final boolean TIMESTAMP_FILES = false;

    private static final String REPORT_BASE_DIR = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";
    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("dd-MMM-yyyy", Locale.ENGLISH);
    private static final String COMMENT_STRING = "This is a test report and is based off COLO829";
    private static final String COMMENT_STRING_CORRECTED = "This is a corrected test report and is based off COLO829";

    private static final String COMMENT_STRING_QC_FAIL = "This is a test QC fail report";
    private static final String COMMENT_STRING_QC_FAIL_CORRECTED = "This is a corrected test QC fail report";

    @Test
    public void canGeneratePatientReportForCOLO829() throws IOException {
        AnalysedPatientReport colo829Report = ExampleAnalysisTestFactory.buildCOLO829(false, COMMENT_STRING);

        CFReportWriter writer = new CFReportWriter(WRITE_TO_PDF);
        writer.writeAnalysedPatientReport(colo829Report, testReportFilePath(colo829Report));
    }

    @Test
    public void canGeneratePatientReportForCOLO829Corrected() throws IOException {
        AnalysedPatientReport colo829Report = ExampleAnalysisTestFactory.buildCOLO829(true, COMMENT_STRING_CORRECTED);

        CFReportWriter writer = new CFReportWriter(WRITE_TO_PDF);
        writer.writeAnalysedPatientReport(colo829Report, testReportFilePath(colo829Report));
    }

    @Test
    public void canGeneratePatientReportForCPCTSample() throws IOException {
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledIn("CPCT01_FULL");

        CFReportWriter writer = new CFReportWriter(WRITE_TO_PDF);
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
    }

    @Test
    public void canGeneratePatientReportForCORESample() throws IOException {
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledIn("CORE01_FULL");

        CFReportWriter writer = new CFReportWriter(WRITE_TO_PDF);
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
    }

    @Test
    public void canGeneratePatientReportForWIDESample() throws IOException {
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledIn("WIDE01_FULL");

        CFReportWriter writer = new CFReportWriter(WRITE_TO_PDF);
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
    }

    @Test
    public void canGeneratePatientReportForBelowDetectionSample() throws IOException {
        AnalysedPatientReport patientReport =
                ExampleAnalysisTestFactory.buildAnalysisWithAllTablesForBelowDetectionLimitSample("CPCT01_NO_TUMOR");

        CFReportWriter writer = new CFReportWriter(WRITE_TO_PDF);
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath(patientReport));
    }

    @Test
    public void canGenerateLowDNAYieldFailReport() throws IOException {
        generateQCFailCPCTReport("CPCT01", "60%", QCFailReason.LOW_DNA_YIELD, false, COMMENT_STRING_QC_FAIL);
    }

    @Test
    public void canGenerateLowDNAYieldFailReportCorrected() throws IOException {
        generateQCFailCPCTReport("CPCT01", "60%", QCFailReason.LOW_DNA_YIELD, true, COMMENT_STRING_QC_FAIL_CORRECTED);
    }

    @Test
    public void canGenerateBelowDetectionThresholdWithoutGenomicAlterationsFailReport() throws IOException {
        generateQCFailCPCTReport("CPCT02", "60%", QCFailReason.BELOW_DETECTION_THRESHOLD, false, COMMENT_STRING_QC_FAIL);
    }

    @Test
    public void canGenerateLabFailureFailReport() throws IOException {
        generateQCFailCPCTReport("CPCT03", "60%", QCFailReason.LAB_FAILURE, false, COMMENT_STRING_QC_FAIL);
    }

    @Test
    public void canGenerateInsufficientTissueFailReport() throws IOException {
        generateQCFailCPCTReport("CPCT04", "60%", QCFailReason.INSUFFICIENT_TISSUE, false, COMMENT_STRING_QC_FAIL);
    }

    @Test
    public void canGeneratePostAnalysisFailReport() throws IOException {
        generateQCFailCPCTReport("CPCT05", "60%", QCFailReason.POST_ANALYSIS_FAIL, false, COMMENT_STRING_QC_FAIL);
    }

    @Test
    public void canGenerateShallowSeqLowTumorPurityFailReport() throws IOException {
        generateQCFailCPCTReport("CPCT06", "15%", QCFailReason.SHALLOW_SEQ_LOW_PURITY, false, COMMENT_STRING_QC_FAIL);
    }

    @Test
    public void canGenerateShallowSeqLowTumorPurityFailReportCORE() throws IOException {
        generateQCFailCPCTReport("CORE01", "15%", QCFailReason.SHALLOW_SEQ_LOW_PURITY, false, COMMENT_STRING_QC_FAIL);
    }

    @Test
    public void canGenerateShallowSeqLowTumorPurityFailReportWIDE() throws IOException {
        generateQCFailCPCTReport("WIDE01", "15%", QCFailReason.SHALLOW_SEQ_LOW_PURITY, false, COMMENT_STRING_QC_FAIL);
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

    private static void generateQCFailCPCTReport(@NotNull String sampleId, @NotNull String shallowSeqPurity, @NotNull QCFailReason reason,
            boolean correctedReport, @NotNull String comments) throws IOException {
        SampleMetadata sampleMetadata = ImmutableSampleMetadata.builder()
                .refSampleId("x")
                .refSampleBarcode("FR12123488")
                .tumorSampleId(sampleId)
                .tumorSampleBarcode("FR12345678")
                .build();

        SampleReport sampleReport = ImmutableSampleReport.builder()
                .sampleMetadata(sampleMetadata)
                .patientTumorLocation(ImmutablePatientTumorLocation.of(Strings.EMPTY, "Skin", "Melanoma"))
                .refArrivalDate(LocalDate.parse("10-Jan-2020", DATE_FORMATTER))
                .tumorArrivalDate(LocalDate.parse("05-Jan-2020", DATE_FORMATTER))
                .shallowSeqPurityString(shallowSeqPurity)
                .labProcedures("PREP013V23-QC037V20-SEQ008V25")
                .cohort("TEST")
                .projectName("TEST-001-002")
                .submissionId("SUBM")
                .hospitalContactData(createTestHospitalContactData())
                .hospitalPatientId("HOSP1")
                .hospitalPathologySampleId("PA1")
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

        String filename = testReportFilePath(patientReport);

        CFReportWriter writer = new CFReportWriter(WRITE_TO_PDF);
        writer.writeQCFailReport(patientReport, filename);
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
