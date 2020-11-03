package com.hartwig.hmftools.patientreporter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import com.google.gson.Gson;
import com.hartwig.hmftools.common.clinical.PatientTumorLocation;
import com.hartwig.hmftools.common.clinical.PatientTumorLocationFile;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
import com.hartwig.hmftools.patientreporter.cfreport.CFReportWriter;
import com.hartwig.hmftools.patientreporter.qcfail.ImmutableQCFailReportData;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReport;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReportData;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReporter;
import com.hartwig.hmftools.patientreporter.reportingdb.ReportingDb;
import com.hartwig.hmftools.protect.variants.germline.GermlineReportingFile;
import com.hartwig.hmftools.protect.variants.germline.GermlineReportingModel;

import org.apache.commons.cli.CommandLine;
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

    // Uncomment this line when generating an example report using PDFWriterTest
    //                public static final String VERSION = "7.17";

    public static void main(@NotNull String[] args) throws IOException {
        Options options = PatientReporterConfig.createOptions();

        PatientReporterConfig config = null;
        try {
            config = PatientReporterConfig.createConfig(createCommandLine(options, args));
        } catch (ParseException exception) {
            LOGGER.warn(exception);
            new HelpFormatter().printHelp("PatientReporter", options);
            System.exit(1);
        }

        LOGGER.info("Running patient reporter v{}", VERSION);
        SampleMetadata sampleMetadata = buildSampleMetadata(config);
        GermlineReportingModel germlineReportingModel = GermlineReportingFile.buildFromTsv(config.germlineReportingTsv());

        ReportWriter reportWriter = CFReportWriter.createProductionReportWriter(germlineReportingModel);
        if (config.qcFail()) {
            LOGGER.info("Generating qc-fail report");
            QCFailReporter reporter = new QCFailReporter(buildBaseReportData(config));
            QCFailReport report = reporter.run(config.qcFailReason(),
                    sampleMetadata,
                    config.purplePurityTsv(),
                    config.purpleQcFile(),
                    config.comments(),
                    config.correctedReport());

            String outputFilePath = generateOutputFilePathForPatientReport(config.outputDirReport(), report);
            reportWriter.writeQCFailReport(report, outputFilePath);

            if (!config.onlyCreatePDF()) {
                LOGGER.debug("Updating additional files and databases");

                writeReportDataToJson(config.outputDirData(),
                        report.sampleReport().tumorSampleId(),
                        report.sampleReport().tumorSampleBarcode(),
                        report);

                ReportingDb.addQCFailReportToReportingDb(config.reportingDbTsv(), report);
            }
        } else {
            LOGGER.info("Generating patient report");
            AnalysedPatientReporter reporter = new AnalysedPatientReporter(buildAnalysedReportData(config));

            AnalysedPatientReport report = reporter.run(sampleMetadata,
                    config.purplePurityTsv(),
                    config.purpleQcFile(),
                    config.purpleDriverCatalogTsv(),
                    config.purpleSomaticVariantVcf(),
                    config.bachelorTsv(),
                    config.linxFusionTsv(),
                    config.linxBreakendTsv(),
                    config.linxViralInsertionTsv(),
                    config.linxDriversTsv(),
                    config.chordPredictionTxt(),
                    config.circosFile(),
                    config.comments(),
                    config.correctedReport());

            String outputFilePath = generateOutputFilePathForPatientReport(config.outputDirReport(), report);
            reportWriter.writeAnalysedPatientReport(report, outputFilePath);

            if (!config.onlyCreatePDF()) {
                LOGGER.debug("Updating additional files and databases");

                writeReportDataToJson(config.outputDirData(),
                        report.sampleReport().sampleMetadata().tumorSampleId(),
                        report.sampleReport().sampleMetadata().tumorSampleBarcode(),
                        report);

                ReportingDb.addAnalysedReportToReportingDb(config.reportingDbTsv(), report);
            }
        }
    }

    private static void writeReportDataToJson(@NotNull String outputDirData, @NotNull String tumorSampleId, @NotNull String tumorBarcode,
            @NotNull PatientReport report) throws IOException {
        String outputFileData = outputDirData + File.separator + tumorSampleId + "_" + tumorBarcode + ".json";
        Gson gson = new Gson();
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFileData));
        writer.write(gson.toJson(report));
        writer.close();
        LOGGER.info("Created json file at {} ", outputFileData);
    }

    @NotNull
    private static String generateOutputFilePathForPatientReport(@NotNull String reportDirectory, @NotNull PatientReport patientReport) {
        return reportDirectory + File.separator + OutputFileUtil.generateOutputFileNameForReport(patientReport);
    }

    @NotNull
    private static SampleMetadata buildSampleMetadata(@NotNull PatientReporterConfig config) {
        SampleMetadata sampleMetadata = ImmutableSampleMetadata.builder()
                .refSampleId(config.refSampleId())
                .refSampleBarcode(config.refSampleBarcode())
                .tumorSampleId(config.tumorSampleId())
                .tumorSampleBarcode(config.tumorSampleBarcode())
                .build();

        LOGGER.info("Printing sample meta data for {}", sampleMetadata.tumorSampleId());
        LOGGER.info(" Tumor sample barcode: {}", sampleMetadata.tumorSampleBarcode());
        LOGGER.info(" Ref sample: {}", sampleMetadata.refSampleId());
        LOGGER.info(" Ref sample barcode: {}", sampleMetadata.refSampleBarcode());

        return sampleMetadata;
    }

    @NotNull
    private static QCFailReportData buildBaseReportData(@NotNull PatientReporterConfig config) throws IOException {
        String tumorLocationTsv = config.tumorLocationTsv();

        List<PatientTumorLocation> patientTumorLocations = PatientTumorLocationFile.read(tumorLocationTsv);
        LOGGER.info("Loaded tumor locations for {} patients from {}", patientTumorLocations.size(), tumorLocationTsv);

        String limsDirectory = config.limsDir();
        Lims lims = LimsFactory.fromLimsDirectory(limsDirectory);
        LOGGER.info("Loaded LIMS data for {} samples from {}", lims.sampleBarcodeCount(), limsDirectory);

        return ImmutableQCFailReportData.builder()
                .patientTumorLocations(patientTumorLocations)
                .limsModel(lims)
                .signaturePath(config.signature())
                .logoRVAPath(config.rvaLogo())
                .logoCompanyPath(config.companyLogo())
                .build();
    }

    @NotNull
    private static AnalysedReportData buildAnalysedReportData(@NotNull PatientReporterConfig config) throws IOException {
        return AnalysedReportDataLoader.buildFromFiles(buildBaseReportData(config),
                config.knowledgebaseDir(),
                config.germlineReportingTsv(),
                config.sampleSummaryTsv());
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull Options options, @NotNull String... args) throws ParseException {
        return new DefaultParser().parse(options, args);
    }
}
