package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.hospital.HospitalModel;
import com.hartwig.hmftools.common.hospital.HospitalModelFactory;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
import com.hartwig.hmftools.common.lims.LimsStudy;
import com.hartwig.hmftools.common.lims.LimsWide;
import com.hartwig.hmftools.common.lims.LimsWideFile;
import com.hartwig.hmftools.patientreporter.cfreport.CFReportWriter;
import com.hartwig.hmftools.patientreporter.qcfail.ImmutableQCFailReportData;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReport;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReportData;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReporter;
import com.hartwig.hmftools.patientreporter.reportingdb.ReportingDb;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class PatientReporterApplication {

    private static final Logger LOGGER = LogManager.getLogger(PatientReporterApplication.class);

    public static final String VERSION = PatientReporterApplication.class.getPackage().getImplementationVersion();

    // Uncomment this line when generating an example report using PDFWriterTest
    //            public static final String VERSION = "7.11";

    // General params needed for every report
    private static final String REF_SAMPLE_ID = "ref_sample_id";
    private static final String REF_SAMPLE_BARCODE = "ref_sample_barcode";
    private static final String TUMOR_SAMPLE_ID = "tumor_sample_id";
    private static final String TUMOR_SAMPLE_BARCODE = "tumor_sample_barcode";
    private static final String OUTPUT_DIRECTORY = "output_dir";

    private static final String REPORTING_DB_TSV = "reporting_db_tsv";
    private static final String TUMOR_LOCATION_CSV = "tumor_location_csv";
    private static final String LIMS_DIRECTORY = "lims_dir";
    private static final String HOSPITAL_DIRECTORY = "hospital_dir";
    private static final String CONTACT_WIDE_TSV = "contact_wide_tsv";

    private static final String RVA_LOGO = "rva_logo";
    private static final String COMPANY_LOGO = "company_logo";
    private static final String SIGNATURE = "signature";

    // Params specific for QC Fail reports
    private static final String QC_FAIL = "qc_fail";
    private static final String QC_FAIL_REASON = "qc_fail_reason";

    // Params specific for actual patient reports
    private static final String PURPLE_PURITY_TSV = "purple_purity_tsv";
    private static final String PURPLE_QC_FILE = "purple_qc_file";
    private static final String PURPLE_GENE_CNV_TSV = "purple_gene_cnv_tsv";
    private static final String SOMATIC_VARIANT_VCF = "somatic_variant_vcf";
    private static final String BACHELOR_TSV = "bachelor_tsv";
    private static final String LINX_FUSION_TSV = "linx_fusion_tsv";
    private static final String LINX_DISRUPTION_TSV = "linx_disruption_tsv";
    private static final String LINX_VIRAL_INSERTION_TSV = "linx_viral_insertion_tsv";
    private static final String LINX_DRIVERS_TSV = "linx_drivers_tsv";
    private static final String CHORD_PREDICTION_TXT = "chord_prediction_txt";
    private static final String CIRCOS_FILE = "circos_file";

    private static final String KNOWLEDGEBASE_DIRECTORY = "knowledgebase_dir";
    private static final String GERMLINE_GENES_CSV = "germline_genes_csv";
    private static final String SAMPLE_SUMMARY_TSV = "sample_summary_tsv";

    // Some additional optional params and flags
    private static final String COMMENTS = "comments";
    private static final String CORRECTED_REPORT = "corrected_report";
    private static final String UNOFFICIAL_REPORT = "unofficial_report";
    private static final String LOG_DEBUG = "log_debug";

    public static void main(final String... args) throws ParseException, IOException {
        final Options options = PatientReporterConfig.createOptions();

        CommandLine cmd = createCommandLine(options, args);
        PatientReporterConfig config = PatientReporterConfig.createConfig(VERSION, cmd);


        if (!PatientReporterConfig.validInputForBaseReport(cmd)) {
            printUsageAndExit(options);
        }

        if (config.logDebug()) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        LOGGER.info("Running patient reporter v{}", VERSION);
        SampleMetadata sampleMetadata = buildSampleMetadata(cmd);

        ReportWriter reportWriter = CFReportWriter.createProductionReportWriter();
        if (config.qcFail() && PatientReporterConfig.validInputForQCFailReport(cmd)) {
            LOGGER.info("Generating qc-fail report");
            QCFailReason reason = QCFailReason.fromIdentifier(config.qcFailReason());
            QCFailReporter reporter = new QCFailReporter(buildQCFailReportData(cmd));
            QCFailReport report = reporter.run(sampleMetadata, reason, config.comments(), config.correctedReport());
            String outputFilePath = generateOutputFilePathForPatientReport(config.outputDir(), report);
            reportWriter.writeQCFailReport(report, outputFilePath);

            ReportingDb.addQCFailReportToReportingDb(config.reportingDbTsv(), report);
        } else if (PatientReporterConfig.validInputForAnalysedReport(cmd)) {
            LOGGER.info("Generating patient report");
            AnalysedPatientReporter reporter = new AnalysedPatientReporter(buildAnalysedReportData(cmd));

            AnalysedPatientReport report = reporter.run(sampleMetadata,
                    cmd.getOptionValue(PURPLE_PURITY_TSV),
                    cmd.getOptionValue(PURPLE_QC_FILE),
                    cmd.getOptionValue(PURPLE_GENE_CNV_TSV),
                    cmd.getOptionValue(SOMATIC_VARIANT_VCF),
                    cmd.getOptionValue(BACHELOR_TSV),
                    cmd.getOptionValue(LINX_FUSION_TSV),
                    cmd.getOptionValue(LINX_DISRUPTION_TSV),
                    cmd.getOptionValue(LINX_VIRAL_INSERTION_TSV),
                    cmd.getOptionValue(LINX_DRIVERS_TSV),
                    cmd.getOptionValue(CHORD_PREDICTION_TXT),
                    cmd.getOptionValue(CIRCOS_FILE),
                    cmd.getOptionValue(COMMENTS),
                    cmd.hasOption(CORRECTED_REPORT),
                    cmd.hasOption(UNOFFICIAL_REPORT));
            String outputFilePath = generateOutputFilePathForPatientReport(cmd.getOptionValue(OUTPUT_DIRECTORY), report);
            reportWriter.writeAnalysedPatientReport(report, outputFilePath);

            ReportingDb.addSequenceReportToReportingDb(cmd.getOptionValue(REPORTING_DB_TSV), report);
        } else {
            printUsageAndExit(options);
        }
    }

    @NotNull
    private static String generateOutputFilePathForPatientReport(@NotNull String reportDirectory, @NotNull PatientReport patientReport) {
        SampleReport sampleReport = patientReport.sampleReport();
        LimsStudy study = LimsStudy.fromSampleId(sampleReport.tumorSampleId());

        String filePrefix = study == LimsStudy.CORE
                ? sampleReport.tumorSampleId() + "_" + sampleReport.hospitalPatientId().replace(" ", "_")
                : sampleReport.tumorSampleId();

        String fileSuffix = patientReport.isCorrectedReport() ? "_corrected.pdf" : ".pdf";

        return reportDirectory + File.separator + filePrefix + "_hmf_report" + fileSuffix;
    }

    @NotNull
    private static SampleMetadata buildSampleMetadata(@NotNull CommandLine cmd) {
        SampleMetadata sampleMetadata = ImmutableSampleMetadata.builder()
                .refSampleId(cmd.getOptionValue(REF_SAMPLE_ID))
                .refSampleBarcode(cmd.getOptionValue(REF_SAMPLE_BARCODE))
                .tumorSampleId(cmd.getOptionValue(TUMOR_SAMPLE_ID))
                .tumorSampleBarcode(cmd.getOptionValue(TUMOR_SAMPLE_BARCODE))
                .build();

        LOGGER.info("Printing sample meta data for {}", sampleMetadata.tumorSampleId());
        LOGGER.info(" Tumor sample barcode: {}", sampleMetadata.tumorSampleBarcode());
        LOGGER.info(" Ref sample: {}", sampleMetadata.refSampleId());
        LOGGER.info(" Ref sample barcode: {}", sampleMetadata.refSampleBarcode());

        return sampleMetadata;
    }

    @NotNull
    private static QCFailReportData buildQCFailReportData(@NotNull CommandLine cmd) throws IOException {
        String tumorLocationCsv = cmd.getOptionValue(TUMOR_LOCATION_CSV);
        List<PatientTumorLocation> patientTumorLocations = PatientTumorLocation.readRecords(tumorLocationCsv);
        LOGGER.info("Loaded tumor locations for {} patients from {}", patientTumorLocations.size(), tumorLocationCsv);

        String limsDirectory = cmd.getOptionValue(LIMS_DIRECTORY);
        Lims lims = LimsFactory.fromLimsDirectory(limsDirectory);
        LOGGER.info("Loaded LIMS data for {} samples from {}", lims.sampleBarcodeCount(), limsDirectory);

        String hospitalsDirectory = cmd.getOptionValue(HOSPITAL_DIRECTORY);
        HospitalModel hospitalModel = HospitalModelFactory.fromHospitalDirectory(hospitalsDirectory);
        LOGGER.info("Loaded data for {} hospitals from {}", hospitalModel.hospitalCount(), hospitalsDirectory);

        LOGGER.info("Reading lims wide file {}", cmd.getOptionValue(CONTACT_WIDE_TSV));
        LimsWide limsWide = LimsWideFile.read(cmd.getOptionValue(CONTACT_WIDE_TSV));

        return ImmutableQCFailReportData.builder()
                .patientTumorLocations(patientTumorLocations)
                .limsModel(lims)
                .limsWideModel(limsWide)
                .hospitalModel(hospitalModel)
                .signaturePath(cmd.getOptionValue(SIGNATURE))
                .logoRVAPath(cmd.getOptionValue(RVA_LOGO))
                .logoCompanyPath(cmd.getOptionValue(COMPANY_LOGO))
                .build();
    }

    @NotNull
    private static AnalysedReportData buildAnalysedReportData(@NotNull CommandLine cmd) throws IOException {
        return AnalysedReportDataLoader.buildFromFiles(buildQCFailReportData(cmd),
                cmd.getOptionValue(KNOWLEDGEBASE_DIRECTORY),
                cmd.getOptionValue(GERMLINE_GENES_CSV),
                cmd.getOptionValue(SAMPLE_SUMMARY_TSV));
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull Options options, @NotNull String... args) throws ParseException {
        return new DefaultParser().parse(options, args);
    }

    private static void printUsageAndExit(@NotNull Options options) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Patient-Reporter", options);
        System.exit(1);
    }
}
