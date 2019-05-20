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
import com.hartwig.hmftools.patientreporter.qcfail.ImmutableQCFailReporter;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReporter;
import com.hartwig.hmftools.patientreporter.structural.SvAnalyzer;

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
    //        public static final String VERSION = "6.1";

    private static final String RUN_DIRECTORY = "run_dir";

    private static final String LIMS_DIRECTORY = "lims_dir";
    private static final String HOSPITALS_DIRECTORY = "hospital_dir";
    private static final String REPORT_DIRECTORY = "report_dir";
    private static final String KNOWLEDGEBASE_DIRECTORY = "knowledgebase_dir";

    private static final String QC_FAIL = "qc_fail";
    private static final String QC_FAIL_REASON = "qc_fail_reason";
    private static final String QC_FAIL_SAMPLE = "qc_fail_sample";

    private static final String FUSION_CSV = "fusion_csv";
    private static final String DISRUPTION_CSV = "disruption_csv";

    private static final String FASTA_FILE_LOCATION = "fasta_file_location";
    private static final String HIGH_CONFIDENCE_BED = "high_confidence_bed";
    private static final String TUMOR_LOCATION_CSV = "tumor_location_csv";
    private static final String DRUP_GENES_CSV = "drup_genes_csv";
    private static final String HOTSPOT_TSV = "hotspot_tsv";
    private static final String GERMLINE_GENES_CSV = "germline_genes_csv";
    private static final String SUMMARY_SAMPLES_CSV = "summary_samples_csv";

    private static final String RVA_LOGO = "rva_logo";
    private static final String COMPANY_LOGO = "company_logo";
    private static final String SIGNATURE = "signature";

    private static final String COMMENTS = "comments";
    private static final String LOG_DEBUG = "log_debug";

    public static void main(final String... args) throws ParseException, IOException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        if (!validInputForReportWriter(cmd) || !validInputForBaseReportData(cmd)) {
            printUsageAndExit(options);
        }

        if (cmd.hasOption(LOG_DEBUG)) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        LOGGER.info("Running patient reporter v" + VERSION);
        final ReportWriter reportWriter = CFReportWriter.createProductionReportWriter();

        if (cmd.hasOption(QC_FAIL) && validInputForQCFailReport(cmd)) {
            final String sample = cmd.getOptionValue(QC_FAIL_SAMPLE);
            LOGGER.info("Generating qc-fail report for {}", sample);
            final QCFailReason reason = QCFailReason.fromIdentifier(cmd.getOptionValue(QC_FAIL_REASON));
            final QCFailReporter reporter = ImmutableQCFailReporter.of(buildBaseReportData(cmd));

            final QCFailReport report = reporter.run(sample, reason, cmd.getOptionValue(COMMENTS));
            final String outputFilePath = generateOutputFilePathForPatientReport(cmd.getOptionValue(REPORT_DIRECTORY), report);
            reportWriter.writeQCFailReport(report, outputFilePath);
        } else if (validInputForAnalysedSample(cmd)) {
            LOGGER.info("Generating sequence report...");
            final SequencedReportData reporterData = buildReporterData(cmd);
            final PatientReporter reporter = buildReporter(cmd, reporterData);

            final AnalysedPatientReport report =
                    reporter.run(cmd.getOptionValue(RUN_DIRECTORY), cmd.getOptionValue(COMMENTS));
            final String outputFilePath = generateOutputFilePathForPatientReport(cmd.getOptionValue(REPORT_DIRECTORY), report);
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
    private static BaseReportData buildBaseReportData(@NotNull final CommandLine cmd) throws IOException {
        String tumorLocationCsv = cmd.getOptionValue(TUMOR_LOCATION_CSV);
        LOGGER.info("Loading tumor location CSV from {}.", tumorLocationCsv);
        final List<PatientTumorLocation> patientTumorLocations = PatientTumorLocation.readRecords(tumorLocationCsv);
        LOGGER.info(" Loaded tumor locations for {} patients.", patientTumorLocations.size());

        String limsDirectory = cmd.getOptionValue(LIMS_DIRECTORY);
        LOGGER.info("Loading LIMS database from {}.", limsDirectory);
        final Lims lims = LimsFactory.fromLimsDirectory(limsDirectory);
        LOGGER.info(" Loaded data for {} samples.", lims.sampleCount());

        String hospitalsDirectory = cmd.getOptionValue(HOSPITALS_DIRECTORY);
        LOGGER.info("Loading hospitals from {}.", hospitalsDirectory);
        final HospitalModel hospitalModel = HospitalModelFactory.fromHospitalDirectory(hospitalsDirectory);
        LOGGER.info(" Loaded data for {} hospitals.", hospitalModel.hospitalCount());

        return ImmutableBaseReportData.of(patientTumorLocations,
                lims,
                hospitalModel,
                cmd.getOptionValue(SIGNATURE),
                cmd.getOptionValue(RVA_LOGO),
                cmd.getOptionValue(COMPANY_LOGO));
    }

    @NotNull
    private static SequencedReportData buildReporterData(@NotNull final CommandLine cmd) throws IOException {
        return SequencedReportDataLoader.buildFromFiles(cmd.getOptionValue(KNOWLEDGEBASE_DIRECTORY),
                cmd.getOptionValue(DRUP_GENES_CSV),
                cmd.getOptionValue(HOTSPOT_TSV),
                cmd.getOptionValue(FASTA_FILE_LOCATION),
                cmd.getOptionValue(HIGH_CONFIDENCE_BED),
                cmd.getOptionValue(GERMLINE_GENES_CSV),
                cmd.getOptionValue(SUMMARY_SAMPLES_CSV));
    }

    @NotNull
    private static PatientReporter buildReporter(@NotNull final CommandLine cmd, @NotNull final SequencedReportData sequencedReportData)
            throws IOException {
        final SvAnalyzer svAnalyzer = SvAnalyzer.fromFiles(cmd.getOptionValue(FUSION_CSV), cmd.getOptionValue(DISRUPTION_CSV));

        return new PatientReporter(buildBaseReportData(cmd), sequencedReportData, svAnalyzer);
    }

    private static boolean validInputForAnalysedSample(@NotNull final CommandLine cmd) {
        final String runDirectory = cmd.getOptionValue(RUN_DIRECTORY);
        final String knowledgebaseDirectory = cmd.getOptionValue(KNOWLEDGEBASE_DIRECTORY);
        final String drupGenesCsv = cmd.getOptionValue(DRUP_GENES_CSV);
        final String hotspotTsv = cmd.getOptionValue(HOTSPOT_TSV);
        final String fusionsCsv = cmd.getOptionValue(FUSION_CSV);
        final String disruptionsCsv = cmd.getOptionValue(DISRUPTION_CSV);
        final String germlineGenesCsv = cmd.getOptionValue(GERMLINE_GENES_CSV);
        final String summarySamplesCsv = cmd.getOptionValue(SUMMARY_SAMPLES_CSV);

        final String fastaFileLocation = cmd.getOptionValue(FASTA_FILE_LOCATION);
        final String highConfidenceBed = cmd.getOptionValue(HIGH_CONFIDENCE_BED);

        if (runDirectory == null || !exists(runDirectory) || !isDirectory(runDirectory)) {
            LOGGER.warn(RUN_DIRECTORY + " has to be an existing directory: " + runDirectory);
        } else if (knowledgebaseDirectory == null || !exists(knowledgebaseDirectory) || !isDirectory(knowledgebaseDirectory)) {
            LOGGER.warn(KNOWLEDGEBASE_DIRECTORY + " has to be an existing directory: " + knowledgebaseDirectory);
        } else if (drupGenesCsv == null || !exists(drupGenesCsv)) {
            LOGGER.warn(DRUP_GENES_CSV + " has to be an existing file: " + drupGenesCsv);
        } else if (hotspotTsv == null || !exists(hotspotTsv)) {
            LOGGER.warn(HOTSPOT_TSV + " has to be an existing file: " + hotspotTsv);
        } else if (fastaFileLocation == null || !exists(fastaFileLocation)) {
            LOGGER.warn(FASTA_FILE_LOCATION + " has to be an existing file: " + fastaFileLocation);
        } else if (fusionsCsv == null || !exists(fusionsCsv)) {
            LOGGER.warn(FUSION_CSV + " has to be an existing file: " + fusionsCsv);
        } else if (disruptionsCsv == null || !exists(disruptionsCsv)) {
            LOGGER.warn(DISRUPTION_CSV + " has to be an existing file: " + disruptionsCsv);
        } else if (highConfidenceBed == null || !exists(highConfidenceBed)) {
            LOGGER.warn(HIGH_CONFIDENCE_BED + " has to be an existing file: " + highConfidenceBed);
        } else if (germlineGenesCsv == null || !exists(germlineGenesCsv)) {
            LOGGER.warn(GERMLINE_GENES_CSV + " has to be an existing file: " + germlineGenesCsv);
        } else if (summarySamplesCsv == null || !exists(summarySamplesCsv)) {
            LOGGER.warn(SUMMARY_SAMPLES_CSV + " has to be an existing file: " + summarySamplesCsv);
        } else {
            return true;
        }
        return false;
    }

    private static boolean validInputForQCFailReport(@NotNull final CommandLine cmd) {
        final QCFailReason qcFailReason = QCFailReason.fromIdentifier(cmd.getOptionValue(QC_FAIL_REASON));
        final String qcFailSample = cmd.getOptionValue(QC_FAIL_SAMPLE);
        if (qcFailReason == QCFailReason.UNDEFINED) {
            LOGGER.warn(QC_FAIL_REASON + " has to be low_tumor_percentage, low_dna_yield, post_analysis_fail or shallow_seq.");
        } else if (qcFailSample == null) {
            LOGGER.warn(QC_FAIL_SAMPLE + " has to be provided.");
        } else {
            return true;
        }
        return false;
    }

    private static boolean validInputForReportWriter(@NotNull final CommandLine cmd) {
        final String reportDirectory = cmd.getOptionValue(REPORT_DIRECTORY);

        if (reportDirectory == null || !exists(reportDirectory) || !isDirectory(reportDirectory)) {
            LOGGER.warn(REPORT_DIRECTORY + " has to be an existing directory: " + reportDirectory);
        } else {
            return true;
        }
        return false;
    }

    private static boolean validInputForBaseReportData(@NotNull final CommandLine cmd) {
        final String tumorLocationCsv = cmd.getOptionValue(TUMOR_LOCATION_CSV);
        final String limsDirectory = cmd.getOptionValue(LIMS_DIRECTORY);
        final String hospitalDirectory = cmd.getOptionValue(HOSPITALS_DIRECTORY);

        final String signaturePath = cmd.getOptionValue(SIGNATURE);
        final String rvaLogoPath = cmd.getOptionValue(RVA_LOGO);
        final String companyLogoPath = cmd.getOptionValue(COMPANY_LOGO);

        if (tumorLocationCsv == null || !exists(tumorLocationCsv)) {
            LOGGER.warn(TUMOR_LOCATION_CSV + " has to be an existing file: " + tumorLocationCsv);
        } else if (limsDirectory == null || !exists(limsDirectory) || !isDirectory(limsDirectory)) {
            LOGGER.warn(LIMS_DIRECTORY + " has to be an existing directory: " + limsDirectory);
        } else if (hospitalDirectory == null || !exists(hospitalDirectory) || !isDirectory(hospitalDirectory)) {
            LOGGER.warn(HOSPITALS_DIRECTORY + " has to be an existing directory: " + hospitalDirectory);
        } else if (signaturePath == null || !exists(signaturePath)) {
            LOGGER.warn(SIGNATURE + " has to be an existing file: " + signaturePath);
        } else if (rvaLogoPath == null || !exists(rvaLogoPath)) {
            LOGGER.warn(RVA_LOGO + " has to be an existing file: " + rvaLogoPath);
        } else if (companyLogoPath == null || !exists(companyLogoPath)) {
            LOGGER.warn(COMPANY_LOGO + " has to be an existing file: " + companyLogoPath);
        } else {
            return true;
        }
        return false;
    }

    private static boolean exists(@NotNull final String path) {
        return Files.exists(new File(path).toPath());
    }

    private static boolean isDirectory(@NotNull final String path) {
        return Files.isDirectory(new File(path).toPath());
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(RUN_DIRECTORY, true, "Path towards a single run dir where patient reporter will run on.");
        options.addOption(LIMS_DIRECTORY, true, "Path a directory holding the LIMS data");
        options.addOption(HOSPITALS_DIRECTORY, true, "Path towards a directory containing hospital data.");
        options.addOption(REPORT_DIRECTORY, true, "Path to where the PDF reports have to be saved.");
        options.addOption(KNOWLEDGEBASE_DIRECTORY, true, "Path towards a directory holding knowledgebase output files.");
        options.addOption(QC_FAIL, false, "If set, generates a qc-fail report.");
        options.addOption(QC_FAIL_REASON, true, "Either 'low_tumor_percentage', 'low_dna_yield', 'post_analysis_fail' or 'shallow_seq'");
        options.addOption(QC_FAIL_SAMPLE, true, "In case of qc-fail reports, the name of the sample used.");
        options.addOption(FUSION_CSV, true, "Path towards a CSV of fusions of the sample");
        options.addOption(DISRUPTION_CSV, true, "Path towards a CSV of disruptions of the sample");
        options.addOption(FASTA_FILE_LOCATION, true, "Path towards the FASTA file containing the ref genome.");
        options.addOption(HIGH_CONFIDENCE_BED, true, "Path towards the high confidence BED file.");
        options.addOption(TUMOR_LOCATION_CSV, true, "Path towards the (curated) tumor location csv.");
        options.addOption(DRUP_GENES_CSV, true, "Path towards a CSV containing genes that could potentially indicate inclusion in DRUP.");
        options.addOption(GERMLINE_GENES_CSV, true, "Path towards a CSV containing germline genes which we want to report");
        options.addOption(SUMMARY_SAMPLES_CSV, true, "Path towards a CSV containing the summary of the samples");
        options.addOption(HOTSPOT_TSV, true, "Path towards a TSV containing known hotspot variants.");
        options.addOption(RVA_LOGO, true, "Path towards a image file containing the RVA logo.");
        options.addOption(COMPANY_LOGO, true, "Path towards a image file containing the company logo.");
        options.addOption(SIGNATURE, true, "Path towards a image file containing the signature to be appended at the end of the report.");
        options.addOption(COMMENTS, true, "Additional comments to be added to the report, if any.");
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
