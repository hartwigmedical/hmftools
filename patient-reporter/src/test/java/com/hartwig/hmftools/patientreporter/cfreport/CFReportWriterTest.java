package com.hartwig.hmftools.patientreporter.cfreport;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testReportData;

import java.io.File;
import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.Locale;
import java.util.Optional;

import com.hartwig.hmftools.common.ecrf.projections.ImmutablePatientTumorLocation;
import com.hartwig.hmftools.common.lims.LimsStudy;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.ExampleAnalysisTestFactory;
import com.hartwig.hmftools.patientreporter.ImmutableSampleMetadata;
import com.hartwig.hmftools.patientreporter.ImmutableSampleReport;
import com.hartwig.hmftools.patientreporter.SampleMetadata;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.qcfail.ImmutableQCFailReport;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReport;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailStudy;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class CFReportWriterTest {

    private static final boolean WRITE_TO_PDF = false;
    private static final boolean TIMESTAMP_FILES = false;

    private static final String REPORT_BASE_DIR = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";
    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("dd-MMM-yyyy", Locale.ENGLISH);

    @Test
    public void canGeneratePatientReportForCOLO829() throws IOException {
        AnalysedPatientReport colo829Report = ExampleAnalysisTestFactory.buildCOLO829();

        CFReportWriter writer = new CFReportWriter(WRITE_TO_PDF);
        writer.writeAnalysedPatientReport(colo829Report, testReportFilePath("hmf_colo829_sequence_report.pdf"));
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
                testReportFilePath("hmf_low_dna_yield_report.pdf"));
    }

    @Test
    public void canGenerateInsufficientTissue() throws IOException {
        generateQCFailCPCTReport("CPCT01000001T",
                "60%",
                QCFailReason.INSUFFICIENT_TISSUE,
                testReportFilePath("hmf_insufficient_tissue_report.pdf"));
    }

    @Test
    public void canGeneratePostDNAIsolationFailReport() throws IOException {
        generateQCFailCPCTReport("CPCT01000001T",
                "60%",
                QCFailReason.POST_ANALYSIS_FAIL,
                testReportFilePath("hmf_post_dna_isolation_fail_report.pdf"));
    }

    @Test
    public void canGenerateLowMolecularTumorPercentageCORE() throws IOException {
        generateQCFailCPCTReport("CORE01000001T",
                "15%",
                QCFailReason.SHALLOW_SEQ_LOW_PURITY,
                testReportFilePath("hmf_low_molecular_tumor_percentage_core_report.pdf"));
    }

    @Test
    public void canGenerateLowMolecularTumorPercentageWIDE() throws IOException {
        generateQCFailCPCTReport("WIDE01000001T",
                "15%",
                QCFailReason.SHALLOW_SEQ_LOW_PURITY,
                testReportFilePath("hmf_low_molecular_tumor_percentage_wide_report.pdf"));
    }

    @Test
    public void canGenerateLowMolecularTumorPercentageCPCT() throws IOException {
        generateQCFailCPCTReport("CPCT01000001T",
                "15%",
                QCFailReason.SHALLOW_SEQ_LOW_PURITY,
                testReportFilePath("hmf_low_molecular_tumor_percentage_cpct_report.pdf"));
    }

    private static void generateQCFailCPCTReport(@NotNull String sampleId, @Nullable String shallowSeqPurity, @NotNull QCFailReason reason,
            @NotNull String filename) throws IOException {
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
                .purityShallowSeq(shallowSeqPurity != null ? shallowSeqPurity : "not determined")
                .labProcedures("PREP013V23-QC037V20-SEQ008V25")
                .addressee("HMF Testing Center")
                .hospitalName(Strings.EMPTY)
                .hospitalPIName(Strings.EMPTY)
                .hospitalPIEmail(Strings.EMPTY)
                .cohort("A")
                .projectName("COLO-001-002")
                .requesterName("ContactMe")
                .requesterEmail("contact@me.com")
                .studyRequesterName("contact")
                .studyRequesterEmail("contact@.com")
                .submissionId("ABC")
                .hospitalPatientId("123456")
                .hospitalPathologySampleId("A")
                .build();

        LimsStudy study = LimsStudy.fromSampleId(sampleId);
        QCFailStudy failStudy;
        switch (study) {
            case CORE:
                failStudy = QCFailStudy.CORE;
                break;
            case WIDE:
                failStudy = QCFailStudy.WIDE;
                break;
            case DRUP:
                failStudy = QCFailStudy.DRUP;
                break;
            default:
                failStudy = QCFailStudy.CPCT;
        }
        QCFailReport patientReport = ImmutableQCFailReport.of(sampleReport,
                reason,
                failStudy,
                Optional.empty(),
                false,
                testReportData().signaturePath(),
                testReportData().logoRVAPath(),
                testReportData().logoCompanyPath());

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
