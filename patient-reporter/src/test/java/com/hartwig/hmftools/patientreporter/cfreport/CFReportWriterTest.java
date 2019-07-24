package com.hartwig.hmftools.patientreporter.cfreport;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testReportData;

import java.io.File;
import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.Locale;
import java.util.Optional;

import com.hartwig.hmftools.common.ecrf.projections.ImmutablePatientTumorLocation;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.ExampleAnalysisTestFactory;
import com.hartwig.hmftools.patientreporter.ImmutableSampleReport;
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
        writer.writeAnalysedPatientReport(colo829Report, testReportFilePath("hmf_test_sequence_report.pdf"));
    }

    @Test
    public void canGeneratePatientReportForCompletelyFilledInReport() throws IOException {
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledIn();

        CFReportWriter writer = new CFReportWriter(WRITE_TO_PDF);
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath("hmf_full_test_sequence_report.pdf"));
    }

    @Test
    public void canGeneratePatientReportForBelowDetectionSamples() throws IOException {
        AnalysedPatientReport patientReport = ExampleAnalysisTestFactory.buildAnalysisWithAllTablesForBelowDetectionLimitSample();

        CFReportWriter writer = new CFReportWriter(WRITE_TO_PDF);
        writer.writeAnalysedPatientReport(patientReport, testReportFilePath("hmf_below_detection_limit_sequence_report.pdf"));
    }

    @Test
    public void canGenerateLowTumorPercentageReport() throws IOException {
        generateQCFailCPCTReport("10%", null, QCFailReason.LOW_TUMOR_PERCENTAGE, testReportFilePath("hmf_low_tumor_percentage_report.pdf"));
    }

    @Test
    public void canGenerateLowDNAYieldReport() throws IOException {
        generateQCFailCPCTReport("60%", null, QCFailReason.LOW_DNA_YIELD, testReportFilePath("hmf_low_dna_yield_report.pdf"));
    }

    @Test
    public void canGenerateInsufficientTissue() throws IOException {
        generateQCFailCPCTReport("60%", null, QCFailReason.INSUFFICIENT_TISSUE, testReportFilePath("hmf_insufficient_tissue_report.pdf"));
    }

    @Test
    public void canGeneratePostDNAIsolationFailReport() throws IOException {
        generateQCFailCPCTReport("60%",
                null,
                QCFailReason.POST_ANALYSIS_FAIL,
                testReportFilePath("hmf_post_dna_isolation_fail_report.pdf"));
    }

    @Test
    public void canGenerateLowMolecularTumorPercentage() throws IOException {
        generateQCFailCPCTReport(null,
                "15%",
                QCFailReason.SHALLOW_SEQ_LOW_PURITY,
                testReportFilePath("hmf_low_molecular_tumor_percentage_report.pdf"));
    }

    private static void generateQCFailCPCTReport(@Nullable String pathologyTumorPercentage, @Nullable String shallowSeqPurity,
            @NotNull QCFailReason reason, @NotNull String filename) throws IOException {
        SampleReport sampleReport = ImmutableSampleReport.builder()
                .sampleId("CPCT02991111T")
                .patientTumorLocation(ImmutablePatientTumorLocation.of("CPCT02991111", "Skin", "Melanoma"))
                .refBarcode("FR12123488")
                .tumorBarcode("FR12345678")
                .refArrivalDate(LocalDate.parse("10-Jan-2019", DATE_FORMATTER))
                .tumorArrivalDate(LocalDate.parse("05-Jan-2019", DATE_FORMATTER))
                .purityShallowSeq(shallowSeqPurity != null ? shallowSeqPurity : "not determined")
                .pathologyTumorPercentage(pathologyTumorPercentage != null ? pathologyTumorPercentage : "not determined")
                .labProcedures("PREP013V23-QC037V20-SEQ008V25")
                .addressee("HMF Testing Center")
                .hospitalName(Strings.EMPTY)
                .hospitalPIName(Strings.EMPTY)
                .hospitalPIEmail(Strings.EMPTY)
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
                Optional.empty(),
                testReportData().signaturePath(),
                testReportData().logoRVAPath(),
                testReportData().logoCompanyPath());

        CFReportWriter writer = new CFReportWriter(WRITE_TO_PDF);
        writer.writeQCFailReport(patientReport, filename);
    }

    @NotNull
    private static String testReportFilePath(@NotNull String filename) {
        if (TIMESTAMP_FILES) {
            int extensionStart = filename.lastIndexOf('.');
            filename = filename.substring(0, extensionStart) + "_" + System.currentTimeMillis() + filename.substring(extensionStart);
        }
        return REPORT_BASE_DIR + File.separator + filename;
    }
}
