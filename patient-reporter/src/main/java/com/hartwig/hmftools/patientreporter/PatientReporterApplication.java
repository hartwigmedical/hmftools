package com.hartwig.hmftools.patientreporter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.gson.GsonBuilder;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumor;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumorFile;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
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
    //                public static final String VERSION = "7.23";

    public static void main(@NotNull String[] args) throws IOException {
        LOGGER.info("Running patient reporter v{}", VERSION);

        Options options = PatientReporterConfig.createOptions();

        PatientReporterConfig config = null;
        try {
            config = PatientReporterConfig.createConfig(new DefaultParser().parse(options, args));
        } catch (ParseException exception) {
            LOGGER.warn(exception);
            new HelpFormatter().printHelp("PatientReporter", options);
            System.exit(1);
        }

        new PatientReporterApplication(config).run();
    }

    @NotNull
    private final PatientReporterConfig config;

    private PatientReporterApplication(@NotNull final PatientReporterConfig config) {
        this.config = config;
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
        AnalysedPatientReporter reporter = new AnalysedPatientReporter(reportData);

        AnalysedPatientReport report = reporter.run(sampleMetadata, config);

        ReportWriter reportWriter = CFReportWriter.createProductionReportWriter();

        String outputFilePath = generateOutputFilePathForPatientReport(config.outputDirReport(), report);
        reportWriter.writeAnalysedPatientReport(report, outputFilePath);

        if (!config.onlyCreatePDF()) {
            LOGGER.debug("Updating reporting db and writing report data");

            writeReportDataToJson(report);

            new ReportingDb(config.reportingDbTsv()).appendAnalysedReport(report);
        }
    }

    private void generateQCFail(@NotNull SampleMetadata sampleMetadata) throws IOException {
        QCFailReporter reporter = new QCFailReporter(buildBaseReportData(config));
        QCFailReport report = reporter.run(config.qcFailReason(),
                sampleMetadata,
                config.purplePurityTsv(),
                config.purpleQcFile(),
                config.comments(),
                config.isCorrectedReport());
        LOGGER.info("Cohort of this sample is: {}", report.sampleReport().cohort().cohortId());

        ReportWriter reportWriter = CFReportWriter.createProductionReportWriter();
        String outputFilePath = generateOutputFilePathForPatientReport(config.outputDirReport(), report);
        reportWriter.writeQCFailReport(report, outputFilePath);

        if (!config.onlyCreatePDF()) {
            LOGGER.debug("Updating reporting db and writing report data");

            writeReportDataToJson(report);

            new ReportingDb(config.reportingDbTsv()).appendQCFailReport(report);
        }
    }

    private void writeReportDataToJson(@NotNull PatientReport report) throws IOException {
        String outputFileData = config.outputDirData() + File.separator + OutputFileUtil.generateOutputFileNameForJson(report);
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFileData));
        writer.write(convertToJson(report));
        writer.close();
        LOGGER.info(" Created report data json file at {} ", outputFileData);
    }

    @VisibleForTesting
    @NotNull
    static String convertToJson(@NotNull PatientReport report) {
        return new GsonBuilder().serializeNulls()
                .serializeSpecialFloatingPointValues()
                .setPrettyPrinting()
                .disableHtmlEscaping()
                .create()
                .toJson(report);
    }

    @NotNull
    private static String generateOutputFilePathForPatientReport(@NotNull String outputDirReport, @NotNull PatientReport patientReport) {
        return outputDirReport + File.separator + OutputFileUtil.generateOutputFileNameForPdfReport(patientReport);
    }

    @NotNull
    private static SampleMetadata buildSampleMetadata(@NotNull PatientReporterConfig config) {
        SampleMetadata sampleMetadata = ImmutableSampleMetadata.builder()
                .refSampleId(config.refSampleId())
                .refSampleRunId(config.refSampleRunId())
                .refSampleBarcode(config.refSampleBarcode())
                .tumorSampleId(config.tumorSampleId())
                .tumorSampleRunId(config.tumorSampleRunId())
                .tumorSampleBarcode(config.tumorSampleBarcode())
                .build();

        LOGGER.info("Printing sample meta data for {}", sampleMetadata.tumorSampleId());
        LOGGER.info(" Tumor sample run id: {}", sampleMetadata.tumorSampleRunId());
        LOGGER.info(" Tumor sample barcode: {}", sampleMetadata.tumorSampleBarcode());
        LOGGER.info(" Ref sample: {}", sampleMetadata.refSampleId());
        LOGGER.info(" Ref sample run id: {}", sampleMetadata.refSampleRunId());
        LOGGER.info(" Ref sample barcode: {}", sampleMetadata.refSampleBarcode());

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
                .build();
    }

    @NotNull
    private static AnalysedReportData buildAnalysedReportData(@NotNull PatientReporterConfig config) throws IOException {
        return AnalysedReportDataLoader.buildFromFiles(buildBaseReportData(config),
                config.germlineReportingTsv(),
                config.sampleSummaryTsv(),
                config.taxonomyDbTsv(),
                config.virusInterpretationTsv(),
                config.virusBlacklistTsv());
    }
}
