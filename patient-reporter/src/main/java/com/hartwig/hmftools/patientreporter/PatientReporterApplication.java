package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.io.IOException;
import java.time.LocalDate;
import java.util.List;

import com.hartwig.hmftools.common.clinical.PatientPrimaryTumor;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumorFile;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
import com.hartwig.hmftools.common.utils.DataUtil;
import com.hartwig.hmftools.patientreporter.algo.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.algo.AnalysedPatientReporter;
import com.hartwig.hmftools.patientreporter.algo.AnalysedReportData;
import com.hartwig.hmftools.patientreporter.algo.AnalysedReportDataLoader;
import com.hartwig.hmftools.patientreporter.cfreport.CFReportWriter;
import com.hartwig.hmftools.patientreporter.qcfail.ImmutableQCFailReportData;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReport;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReportData;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReporter;
import com.hartwig.hmftools.patientreporter.reportingdb.ReportingDb;

import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class PatientReporterApplication {

    private static final Logger LOGGER = LogManager.getLogger(PatientReporterApplication.class);

    public static final String VERSION = PatientReporterApplication.class.getPackage().getImplementationVersion();

    // Uncomment this line when generating an example report using CFReportWriterTest
    //                public static final String VERSION = "7.25";

    public static void main(@NotNull String[] args) throws IOException {
        LOGGER.info("Running patient reporter v{}", VERSION);

        Options options = PatientReporterConfig.createOptions();

        PatientReporterConfig config = null;
        try {
            config = PatientReporterConfig.createConfig(new DefaultParser().parse(options, args));
        } catch (ParseException exception) {
            LOGGER.warn(exception);
            new HelpFormatter().printHelp("PatientReporter", options);
            throw new IllegalArgumentException("Unexpected error, check inputs");
        }

        new PatientReporterApplication(config, DataUtil.formatDate(LocalDate.now())).run();
    }

    @NotNull
    private final PatientReporterConfig config;
    @NotNull
    private final String reportDate;

    private PatientReporterApplication(@NotNull final PatientReporterConfig config, @NotNull final String reportDate) {
        this.config = config;
        this.reportDate = reportDate;
    }

    private void run() throws IOException {
        SampleMetadata sampleMetadata = buildSampleMetadata(config);

        if (config.qcFail()) {
            LOGGER.info("Generating qc-fail report");
            generateQCFail(sampleMetadata);
        } else {
            LOGGER.info("Generating patient report");
            generateAnalysedReport(sampleMetadata);
        }
    }

    private void generateAnalysedReport(@NotNull SampleMetadata sampleMetadata) throws IOException {
        AnalysedReportData reportData = buildAnalysedReportData(config);
        AnalysedPatientReporter reporter = new AnalysedPatientReporter(reportData, reportDate);

        AnalysedPatientReport report = reporter.run(sampleMetadata, config);

        ReportWriter reportWriter = CFReportWriter.createProductionReportWriter();

        String outputFilePath = generateOutputFilePathForPatientReport(config.outputDirReport(), report);
        reportWriter.writeAnalysedPatientReport(report, outputFilePath);

        if (!config.onlyCreatePDF()) {
            LOGGER.debug("Updating reporting db and writing report data");

            reportWriter.writeJsonAnalysedFile(report, config.outputDirData());

            new ReportingDb().appendAnalysedReport(report, config.outputDirData());
        }
    }

    private void generateQCFail(@NotNull SampleMetadata sampleMetadata) throws IOException {
        QCFailReporter reporter = new QCFailReporter(buildBaseReportData(config), reportDate);
        QCFailReport report = reporter.run(sampleMetadata, config);
        LOGGER.info("Cohort of this sample is: {}", report.sampleReport().cohort().cohortId());

        ReportWriter reportWriter = CFReportWriter.createProductionReportWriter();
        String outputFilePath = generateOutputFilePathForPatientReport(config.outputDirReport(), report);

        reportWriter.writeQCFailReport(report, outputFilePath);

        if (!config.onlyCreatePDF()) {
            LOGGER.debug("Updating reporting db and writing report data");

            reportWriter.writeJsonFailedFile(report, config.outputDirData());

            new ReportingDb().appendQCFailReport(report, config.outputDirReport());
        }
    }

    @NotNull
    private static String generateOutputFilePathForPatientReport(@NotNull String outputDirReport, @NotNull PatientReport patientReport) {
        return outputDirReport + File.separator + OutputFileUtil.generateOutputFileNameForPdfReport(patientReport);
    }

    @NotNull
    private static SampleMetadata buildSampleMetadata(@NotNull PatientReporterConfig config) {
        String sampleNameForReport = config.sampleNameForReport();
        SampleMetadata sampleMetadata = ImmutableSampleMetadata.builder()
                .refSampleId(config.refSampleId())
                .refSampleBarcode(config.refSampleBarcode())
                .tumorSampleId(config.tumorSampleId())
                .tumorSampleBarcode(config.tumorSampleBarcode())
                .sampleNameForReport(sampleNameForReport != null ? sampleNameForReport : config.tumorSampleId())
                .build();

        LOGGER.info("Printing sample meta data for {}", sampleMetadata.tumorSampleId());
        LOGGER.info(" Tumor sample barcode: {}", sampleMetadata.tumorSampleBarcode());
        LOGGER.info(" Ref sample: {}", sampleMetadata.refSampleId());
        LOGGER.info(" Ref sample barcode: {}", sampleMetadata.refSampleBarcode());
        LOGGER.info(" Sample name for report: {}", sampleMetadata.sampleNameForReport());

        return sampleMetadata;
    }

    @NotNull
    private static QCFailReportData buildBaseReportData(@NotNull PatientReporterConfig config) throws IOException {
        String primaryTumorTsv = config.primaryTumorTsv();

        List<PatientPrimaryTumor> patientPrimaryTumors = PatientPrimaryTumorFile.read(primaryTumorTsv);
        LOGGER.info("Loaded primary tumors for {} patients from {}", patientPrimaryTumors.size(), primaryTumorTsv);

        String limsDirectory = config.limsDir();
        Lims lims = LimsFactory.fromLimsDirectory(limsDirectory);
        LOGGER.info("Loaded LIMS data for {} samples from {}", lims.sampleBarcodeCount(), limsDirectory);

        return ImmutableQCFailReportData.builder()
                .patientPrimaryTumors(patientPrimaryTumors)
                .limsModel(lims)
                .signaturePath(config.signature())
                .logoRVAPath(config.rvaLogo())
                .logoCompanyPath(config.companyLogo())
                .udiDi(config.udiDi())
                .build();
    }

    @NotNull
    private static AnalysedReportData buildAnalysedReportData(@NotNull PatientReporterConfig config) throws IOException {
        return AnalysedReportDataLoader.buildFromFiles(buildBaseReportData(config),
                config.germlineReportingTsv(),
                config.sampleSummaryTsv(),
                config.sampleSpecialRemarkTsv());
    }
}