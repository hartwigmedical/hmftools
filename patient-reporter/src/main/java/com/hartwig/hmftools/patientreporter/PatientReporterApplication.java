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
import com.hartwig.hmftools.common.lims.LimsSampleType;
import com.hartwig.hmftools.patientreporter.cfreport.CFReportWriter;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReporter;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
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
    //        public static final String VERSION = "7.0";

    // General params needed for every report
    private static final String TUMOR_SAMPLE = "tumor_sample";
    private static final String OUTPUT_DIRECTORY = "output_dir";

    private static final String LIMS_DIRECTORY = "lims_dir";
    private static final String HOSPITAL_DIRECTORY = "hospital_dir";
    private static final String TUMOR_LOCATION_CSV = "tumor_location_csv";

    private static final String RVA_LOGO = "rva_logo";
    private static final String COMPANY_LOGO = "company_logo";
    private static final String SIGNATURE = "signature";

    // Params specific for QC Fail reports
    private static final String QC_FAIL = "qc_fail";
    private static final String QC_FAIL_REASON = "qc_fail_reason";

    // Params specific for actual patient reports
    private static final String PURPLE_PURITY_TSV = "purple_purity_tsv";
    private static final String PURPLE_GENE_CNV_TSV = "purple_gene_cnv_tsv";
    private static final String SOMATIC_VARIANT_VCF = "somatic_variant_vcf";
    private static final String LINX_FUSION_TSV = "linx_fusion_tsv";
    private static final String LINX_DISRUPTION_TSV = "linx_disruption_tsv";
    private static final String BACHELOR_CSV = "bachelor_csv";
    private static final String CHORD_PREDICTION_FILE = "chord_prediction_file";
    private static final String CIRCOS_FILE = "circos_file";

    private static final String REF_SAMPLE = "ref_sample";
    private static final String KNOWLEDGEBASE_DIRECTORY = "knowledgebase_dir";
    private static final String DRUP_GENES_CSV = "drup_genes_csv";
    private static final String GERMLINE_GENES_CSV = "germline_genes_csv";
    private static final String SAMPLE_SUMMARY_CSV = "sample_summary_csv";
    private static final String FASTA_FILE_LOCATION = "fasta_file_location";
    private static final String HIGH_CONFIDENCE_BED = "high_confidence_bed";

    // Some additional optional params
    private static final String COMMENTS = "comments";
    private static final String LOG_DEBUG = "log_debug";

    public static void main(final String... args) throws ParseException, IOException {
        Options options = createOptions();
        CommandLine cmd = createCommandLine(options, args);

        if (!validInputForReportWriter(cmd) || !validInputForBaseReportData(cmd)) {
            printUsageAndExit(options);
        }

        if (cmd.hasOption(LOG_DEBUG)) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        LOGGER.info("Running patient reporter v" + VERSION);
        ReportWriter reportWriter = CFReportWriter.createProductionReportWriter();

        if (cmd.hasOption(QC_FAIL) && validInputForQCFailReport(cmd)) {
            String tumorSample = cmd.getOptionValue(TUMOR_SAMPLE);
            LOGGER.info("Generating qc-fail report for {}", tumorSample);
            QCFailReason reason = QCFailReason.fromIdentifier(cmd.getOptionValue(QC_FAIL_REASON));
            QCFailReporter reporter = new QCFailReporter(buildBaseReportData(cmd));

            QCFailReport report = reporter.run(tumorSample, reason, cmd.getOptionValue(COMMENTS));
            String outputFilePath = generateOutputFilePathForPatientReport(cmd.getOptionValue(OUTPUT_DIRECTORY), report);
            reportWriter.writeQCFailReport(report, outputFilePath);
        } else if (validInputForAnalysedSample(cmd)) {
            String tumorSample = cmd.getOptionValue(TUMOR_SAMPLE);
            LOGGER.info("Generating patient report for {}", tumorSample);
            SequencedReportData sequencedReportData = buildSequencedReportData(cmd);
            PatientReporter reporter = new PatientReporter(sequencedReportData);

            AnalysedPatientReport report = reporter.run(tumorSample,
                    cmd.getOptionValue(REF_SAMPLE),
                    cmd.getOptionValue(PURPLE_PURITY_TSV),
                    cmd.getOptionValue(PURPLE_GENE_CNV_TSV),
                    cmd.getOptionValue(SOMATIC_VARIANT_VCF),
                    cmd.getOptionValue(LINX_FUSION_TSV),
                    cmd.getOptionValue(LINX_DISRUPTION_TSV),
                    cmd.getOptionValue(BACHELOR_CSV),
                    cmd.getOptionValue(CHORD_PREDICTION_FILE),
                    cmd.getOptionValue(CIRCOS_FILE),
                    cmd.getOptionValue(COMMENTS));
            String outputFilePath = generateOutputFilePathForPatientReport(cmd.getOptionValue(OUTPUT_DIRECTORY), report);
            reportWriter.writeAnalysedPatientReport(report, outputFilePath);
        } else {
            printUsageAndExit(options);
        }
    }

    @NotNull
    private static String generateOutputFilePathForPatientReport(@NotNull String reportDirectory, @NotNull PatientReport patientReport) {
        SampleReport sampleReport = patientReport.sampleReport();
        LimsSampleType type = LimsSampleType.fromSampleId(sampleReport.sampleId());

        String filePrefix =
                type == LimsSampleType.CORE ? sampleReport.sampleId() + "_" + sampleReport.hospitalPatientId() : sampleReport.sampleId();
        return reportDirectory + File.separator + filePrefix + "_hmf_report.pdf";
    }

    @NotNull
    private static BaseReportData buildBaseReportData(@NotNull CommandLine cmd) throws IOException {
        String tumorLocationCsv = cmd.getOptionValue(TUMOR_LOCATION_CSV);
        List<PatientTumorLocation> patientTumorLocations = PatientTumorLocation.readRecords(tumorLocationCsv);
        LOGGER.info("Loaded tumor locations for {} patients from {}", patientTumorLocations.size(), tumorLocationCsv);

        String limsDirectory = cmd.getOptionValue(LIMS_DIRECTORY);
        Lims lims = LimsFactory.fromLimsDirectory(limsDirectory);
        LOGGER.info("Loaded LIMS data for {} samples from {}", lims.sampleCount(), limsDirectory);

        String hospitalsDirectory = cmd.getOptionValue(HOSPITAL_DIRECTORY);
        HospitalModel hospitalModel = HospitalModelFactory.fromHospitalDirectory(hospitalsDirectory);
        LOGGER.info("Loaded data for {} hospitals from {}", hospitalModel.hospitalCount(), hospitalsDirectory);

        return ImmutableQCFailReportData.builder()
                .patientTumorLocations(patientTumorLocations)
                .limsModel(lims)
                .hospitalModel(hospitalModel)
                .signaturePath(cmd.getOptionValue(SIGNATURE))
                .logoRVAPath(cmd.getOptionValue(RVA_LOGO))
                .logoCompanyPath(cmd.getOptionValue(COMPANY_LOGO))
                .build();
    }

    @NotNull
    private static SequencedReportData buildSequencedReportData(@NotNull CommandLine cmd) throws IOException {
        return SequencedReportDataLoader.buildFromFiles(buildBaseReportData(cmd),
                cmd.getOptionValue(KNOWLEDGEBASE_DIRECTORY),
                cmd.getOptionValue(DRUP_GENES_CSV),
                cmd.getOptionValue(FASTA_FILE_LOCATION),
                cmd.getOptionValue(HIGH_CONFIDENCE_BED),
                cmd.getOptionValue(GERMLINE_GENES_CSV),
                cmd.getOptionValue(SAMPLE_SUMMARY_CSV));
    }

    private static boolean validInputForAnalysedSample(@NotNull CommandLine cmd) {
        return fileExists(cmd, PURPLE_PURITY_TSV) && fileExists(cmd, PURPLE_GENE_CNV_TSV) && fileExists(cmd, SOMATIC_VARIANT_VCF)
                && fileExists(cmd, LINX_FUSION_TSV) && fileExists(cmd, LINX_DISRUPTION_TSV) && valueMissingOrFileExists(cmd, BACHELOR_CSV)
                && fileExists(cmd, CHORD_PREDICTION_FILE) && fileExists(cmd, CIRCOS_FILE) && valueExists(cmd, REF_SAMPLE) && dirExists(cmd,
                KNOWLEDGEBASE_DIRECTORY) && fileExists(cmd, DRUP_GENES_CSV) && fileExists(cmd, GERMLINE_GENES_CSV) && fileExists(cmd,
                SAMPLE_SUMMARY_CSV) && fileExists(cmd, FASTA_FILE_LOCATION) && fileExists(cmd, HIGH_CONFIDENCE_BED);
    }

    private static boolean validInputForQCFailReport(@NotNull CommandLine cmd) {
        final QCFailReason qcFailReason = QCFailReason.fromIdentifier(cmd.getOptionValue(QC_FAIL_REASON));
        if (qcFailReason == QCFailReason.UNDEFINED) {
            LOGGER.warn(QC_FAIL_REASON + " has to be low_tumor_percentage, low_dna_yield, post_analysis_fail or shallow_seq_low_purity.");
        } else {
            return true;
        }
        return false;
    }

    private static boolean validInputForReportWriter(@NotNull CommandLine cmd) {
        return valueExists(cmd, TUMOR_SAMPLE) && dirExists(cmd, OUTPUT_DIRECTORY);
    }

    private static boolean validInputForBaseReportData(@NotNull CommandLine cmd) {
        return fileExists(cmd, TUMOR_LOCATION_CSV) && dirExists(cmd, LIMS_DIRECTORY) && dirExists(cmd, HOSPITAL_DIRECTORY)
                && fileExists(cmd, SIGNATURE) && fileExists(cmd, RVA_LOGO) && fileExists(cmd, COMPANY_LOGO);
    }

    private static boolean valueExists(@NotNull CommandLine cmd, @NotNull String param) {
        String value = cmd.getOptionValue(param);
        if (value == null) {
            LOGGER.warn(param + " has to be provided");
            return false;
        }
        return true;
    }

    private static boolean valueMissingOrFileExists(@NotNull CommandLine cmd, @NotNull String param) {
        String value = cmd.getOptionValue(param);

        if (value != null && !Files.exists(new File(value).toPath())) {
            LOGGER.warn(param + " is optional, but when provided it has to be an existing file: " + value);
            return false;
        }

        return true;
    }

    private static boolean fileExists(@NotNull CommandLine cmd, @NotNull String param) {
        String value = cmd.getOptionValue(param);

        if (value == null || !pathExists(value)) {
            LOGGER.warn(param + " has to be an existing file: " + value);
            return false;
        }

        return true;
    }

    private static boolean dirExists(@NotNull CommandLine cmd, @NotNull String param) {
        String value = cmd.getOptionValue(param);

        if (value == null || !pathExists(value) || !pathIsDirectory(value)) {
            LOGGER.warn(param + " has to be an existing directory: " + value);
            return false;
        }

        return true;
    }

    private static boolean pathExists(@NotNull String path) {
        return Files.exists(new File(path).toPath());
    }

    private static boolean pathIsDirectory(@NotNull String path) {
        return Files.isDirectory(new File(path).toPath());
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(TUMOR_SAMPLE, true, "The sample for which a patient report will be generated.");
        options.addOption(OUTPUT_DIRECTORY, true, "Path to where the PDF reports have to be written to.");

        options.addOption(TUMOR_LOCATION_CSV, true, "Path towards the (curated) tumor location CSV.");
        options.addOption(LIMS_DIRECTORY, true, "Path towards the directory holding the LIMS data");
        options.addOption(HOSPITAL_DIRECTORY, true, "Path towards the directory containing hospital data.");

        options.addOption(RVA_LOGO, true, "Path towards a image file containing the RVA logo.");
        options.addOption(COMPANY_LOGO, true, "Path towards a image file containing the company logo.");
        options.addOption(SIGNATURE, true, "Path towards a image file containing the signature to be appended at the end of the report.");

        options.addOption(QC_FAIL, false, "If set, generates a qc-fail report.");
        options.addOption(QC_FAIL_REASON,
                true,
                "Either 'low_tumor_percentage', 'low_dna_yield', 'post_analysis_fail', 'shallow_seq' or 'insufficient_tissue_delivered'");

        options.addOption(PURPLE_PURITY_TSV, true, "Path towards the purple purity TSV.");
        options.addOption(PURPLE_GENE_CNV_TSV, true, "Path towards the purple gene copy number TSV.");
        options.addOption(SOMATIC_VARIANT_VCF, true, "Path towards the somatic variant VCF.");
        options.addOption(LINX_FUSION_TSV, true, "Path towards the linx fusion TSV.");
        options.addOption(LINX_DISRUPTION_TSV, true, "Path towards the linx disruption TSV.");
        options.addOption(BACHELOR_CSV, true, "Path towards the germline germline CSV (optional).");
        options.addOption(CHORD_PREDICTION_FILE, true, "Path towards the CHORD prediction file.");
        options.addOption(CIRCOS_FILE, true, "Path towards the circos file.");

        options.addOption(REF_SAMPLE, true, "The reference sample for the sample for which we are generating a report.");
        options.addOption(KNOWLEDGEBASE_DIRECTORY, true, "Path towards the directory holding knowledgebase output files.");
        options.addOption(FASTA_FILE_LOCATION, true, "Path towards the FASTA file containing the ref genome.");
        options.addOption(HIGH_CONFIDENCE_BED, true, "Path towards the high confidence BED file.");
        options.addOption(DRUP_GENES_CSV, true, "Path towards a CSV containing genes that could potentially indicate inclusion in DRUP.");
        options.addOption(GERMLINE_GENES_CSV, true, "Path towards a CSV containing germline genes which we want to report.");
        options.addOption(SAMPLE_SUMMARY_CSV, true, "Path towards a CSV containing the (clinical) summaries of the samples.");

        options.addOption(COMMENTS, true, "Additional comments to be added to the report (optional).");
        options.addOption(LOG_DEBUG, false, "If provided, set the log level to debug rather than default.");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    private static void printUsageAndExit(@NotNull final Options options) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Patient-Reporter", options);
        System.exit(1);
    }
}
