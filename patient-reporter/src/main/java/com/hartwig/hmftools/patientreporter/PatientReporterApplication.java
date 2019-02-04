package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.hartwig.hmftools.common.center.Center;
import com.hartwig.hmftools.common.center.CenterModel;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
import com.hartwig.hmftools.patientreporter.qcfail.ImmutableNotAnalysableReporter;
import com.hartwig.hmftools.patientreporter.qcfail.NotAnalysableReason;
import com.hartwig.hmftools.patientreporter.qcfail.NotAnalysableReporter;
import com.hartwig.hmftools.patientreporter.qcfail.NotAnalysableStudy;
import com.hartwig.hmftools.patientreporter.report.PDFWriter;
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

import net.sf.dynamicreports.report.exception.DRException;

public class PatientReporterApplication {

    private static final Logger LOGGER = LogManager.getLogger(PatientReporterApplication.class);

    public static final String VERSION = PatientReporterApplication.class.getPackage().getImplementationVersion();

    // Uncomment this line when generating an example report using PDFWriterTest
    //        public static final String VERSION = "5.7";

    private static final String TUMOR_LOCATION_CSV = "tumor_location_csv";
    private static final String LIMS_DIRECTORY = "lims";
    private static final String REPORT_DIRECTORY = "report_dir";
    private static final String RUN_DIRECTORY = "run_dir";
    private static final String NOT_ANALYSABLE = "not_analysable";
    private static final String NOT_ANALYSABLE_REASON = "not_analysable_reason";
    private static final String NOT_ANALYSED_SAMPLE = "not_analysable_sample";

    private static final String DO_REPORT_GERMLINE = "do_report_germline";
    private static final String FUSION_CSV = "fusion_csv";
    private static final String DISRUPTION_CSV = "disruption_csv";

    private static final String KNOWLEDGEBASE_PATH = "knowledgebase_path";
    private static final String DRUP_GENES_CSV = "drup_genes_csv";
    private static final String HOTSPOT_TSV = "hotspot_tsv";
    private static final String LOG_DEBUG = "log_debug";

    private static final String FASTA_FILE_LOCATION = "fasta_file_location";
    private static final String HIGH_CONFIDENCE_BED = "high_confidence_bed";

    private static final String CENTER_CSV = "center_csv";
    private static final String RVA_LOGO = "rva_logo";
    private static final String SIGNATURE = "signature";
    private static final String COMMENTS = "comments";

    public static void main(final String... args) throws ParseException, IOException, DRException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        if (!validInputForReportWriter(cmd) || !validInputForBaseReportData(cmd)) {
            printUsageAndExit(options);
        }

        if (cmd.hasOption(LOG_DEBUG)) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        LOGGER.info("Running patient reporter v" + VERSION);
        final PDFWriter pdfWriter = new PDFWriter(cmd.getOptionValue(REPORT_DIRECTORY));

        if (cmd.hasOption(NOT_ANALYSABLE) && validInputForNonAnalysableReport(cmd)) {
            final String sample = cmd.getOptionValue(NOT_ANALYSED_SAMPLE);
            LOGGER.info("Generating non-analysable report for {}", sample);
            final NotAnalysableReason reason = NotAnalysableReason.fromIdentifier(cmd.getOptionValue(NOT_ANALYSABLE_REASON));
            final NotAnalysableReporter reporter = ImmutableNotAnalysableReporter.of(buildBaseReportData(cmd));

            final NotAnalysedPatientReport report = reporter.run(sample, reason, cmd.getOptionValue(COMMENTS));
            pdfWriter.writeNonSequenceableReport(report);
        } else if (validInputForPatientReporter(cmd)) {
            LOGGER.info("Generating sequence report...");
            final SequencedReportData reporterData = buildReporterData(cmd);
            final PatientReporter reporter = buildReporter(cmd, reporterData);

            final AnalysedPatientReport report =
                    reporter.run(cmd.getOptionValue(RUN_DIRECTORY), cmd.hasOption(DO_REPORT_GERMLINE), cmd.getOptionValue(COMMENTS));
            pdfWriter.writeSequenceReport(report);
        } else {
            printUsageAndExit(options);
        }
    }

    @NotNull
    private static BaseReportData buildBaseReportData(@NotNull final CommandLine cmd) throws IOException {
        String tumorLocationCsv = cmd.getOptionValue(TUMOR_LOCATION_CSV);
        LOGGER.info("Loading ECRF CSV dump from {}.", tumorLocationCsv);
        final List<PatientTumorLocation> patientTumorLocations = PatientTumorLocation.readRecords(tumorLocationCsv);
        LOGGER.info(" Loaded data for {} patients.", patientTumorLocations.size());

        String limsDirectory = cmd.getOptionValue(LIMS_DIRECTORY);
        LOGGER.info("Loading LIMS database from {}.", limsDirectory);
        final Lims lims = LimsFactory.fromLimsDirectory(limsDirectory);
        LOGGER.info(" Loaded data for {} samples.", lims.sampleCount());

        final CenterModel centerModel = Center.readFromCSV(cmd.getOptionValue(CENTER_CSV));
        return ImmutableBaseReportData.of(patientTumorLocations,
                lims,
                centerModel,
                cmd.getOptionValue(SIGNATURE),
                cmd.getOptionValue(RVA_LOGO));
    }

    @NotNull
    private static SequencedReportData buildReporterData(@NotNull final CommandLine cmd) throws IOException {
        return SequencedReportDataLoader.buildFromFiles(cmd.getOptionValue(KNOWLEDGEBASE_PATH),
                cmd.getOptionValue(DRUP_GENES_CSV),
                cmd.getOptionValue(HOTSPOT_TSV),
                cmd.getOptionValue(FASTA_FILE_LOCATION),
                cmd.getOptionValue(HIGH_CONFIDENCE_BED));
    }

    @NotNull
    private static PatientReporter buildReporter(@NotNull final CommandLine cmd, @NotNull final SequencedReportData sequencedReportData)
            throws IOException {
        final SvAnalyzer svAnalyzer =
                SvAnalyzer.fromFiles(cmd.getOptionValue(FUSION_CSV), cmd.getOptionValue(DISRUPTION_CSV));

        return ImmutablePatientReporter.of(buildBaseReportData(cmd), sequencedReportData, svAnalyzer);
    }

    private static boolean validInputForPatientReporter(@NotNull final CommandLine cmd) {
        final String runDirectory = cmd.getOptionValue(RUN_DIRECTORY);
        final String knowledgebasePath = cmd.getOptionValue(KNOWLEDGEBASE_PATH);
        final String drupGenesCsv = cmd.getOptionValue(DRUP_GENES_CSV);
        final String hotspotTsv = cmd.getOptionValue(HOTSPOT_TSV);
        final String fusionsCsv = cmd.getOptionValue(FUSION_CSV);
        final String disruptionsCsv = cmd.getOptionValue(DISRUPTION_CSV);

        final String fastaFileLocation = cmd.getOptionValue(FASTA_FILE_LOCATION);
        final String highConfidenceBed = cmd.getOptionValue(HIGH_CONFIDENCE_BED);

        if (runDirectory == null || !exists(runDirectory) || !isDirectory(runDirectory)) {
            LOGGER.warn(RUN_DIRECTORY + " has to be an existing directory: " + runDirectory);
        } else if (knowledgebasePath == null || !exists(knowledgebasePath) || !isDirectory(knowledgebasePath)) {
            LOGGER.warn(KNOWLEDGEBASE_PATH + " has to be an existing directory: " + knowledgebasePath);
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
        } else {
            return true;
        }
        return false;
    }

    private static boolean validInputForNonAnalysableReport(@NotNull final CommandLine cmd) {
        final NotAnalysableReason notAnalysableReason = NotAnalysableReason.fromIdentifier(cmd.getOptionValue(NOT_ANALYSABLE_REASON));
        final String notAnalysedSample = cmd.getOptionValue(NOT_ANALYSED_SAMPLE);
        LOGGER.info("core: " + notAnalysedSample);
        if (notAnalysableReason == NotAnalysableReason.UNDEFINED) {
            LOGGER.warn(NOT_ANALYSABLE_REASON + " has to be low_tumor_percentage, low_dna_yield or post_analysis_fail.");
        } else if (notAnalysedSample == null) {
            LOGGER.warn(NOT_ANALYSED_SAMPLE + " has to be provided.");
        } else if (NotAnalysableStudy.fromSample(notAnalysedSample) == null &&  !notAnalysedSample.contains("CORE")) {
            LOGGER.warn("Could not determine study for sample " + notAnalysedSample);
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
        final String centerCsv = cmd.getOptionValue(CENTER_CSV);
        final String signaturePath = cmd.getOptionValue(SIGNATURE);
        final String rvaLogoPath = cmd.getOptionValue(RVA_LOGO);

        if (tumorLocationCsv == null || !exists(tumorLocationCsv)) {
            LOGGER.warn(TUMOR_LOCATION_CSV + " has to be an existing file: " + tumorLocationCsv);
        } else if (limsDirectory == null || !exists(limsDirectory) || !isDirectory(limsDirectory)) {
            LOGGER.warn(LIMS_DIRECTORY + " has to be an existing directory: " + limsDirectory);
        } else if (centerCsv == null || !exists(centerCsv)) {
            LOGGER.warn(CENTER_CSV + " has to be an existing file: " + centerCsv);
        } else if (signaturePath == null || !exists(signaturePath)) {
            LOGGER.warn(SIGNATURE + " has to be an existing file: " + signaturePath);
        } else if (rvaLogoPath == null || !exists(rvaLogoPath)) {
            LOGGER.warn(RVA_LOGO + " has to be an existing file: " + rvaLogoPath);
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
        options.addOption(TUMOR_LOCATION_CSV, true, "Complete path towards the (curated) tumor location csv.");
        options.addOption(LIMS_DIRECTORY, true, "Complete path a directory holding the LIMS data");
        options.addOption(REPORT_DIRECTORY, true, "Complete path to where the PDF reports have to be saved.");
        options.addOption(RUN_DIRECTORY, true, "Complete path towards a single run dir where patient reporter will run on.");
        options.addOption(NOT_ANALYSABLE, false, "If set, generates a non-analysable report.");
        options.addOption(NOT_ANALYSABLE_REASON, true, "Either 'low_tumor_percentage', 'low_dna_yield' or 'post_analysis_fail'");
        options.addOption(NOT_ANALYSED_SAMPLE, true, "In case of non-sequenceable reports, the name of the sample used.");
        options.addOption(DO_REPORT_GERMLINE, false, "If provided, report germline. Otherwise do not report germline.");
        options.addOption(KNOWLEDGEBASE_PATH, true, "Path towards a directory holding knowledgebase output files.");
        options.addOption(DRUP_GENES_CSV, true, "Path towards a CSV containing genes that could potentially indicate inclusion in DRUP.");
        options.addOption(HOTSPOT_TSV, true, "Path towards a TSV containing known hotspot variants.");
        options.addOption(LOG_DEBUG, false, "If provided, set the log level to debug rather than default.");
        options.addOption(FASTA_FILE_LOCATION, true, "Path towards the FASTA file containing the ref genome.");
        options.addOption(HIGH_CONFIDENCE_BED, true, "Path towards the high confidence BED file.");
        options.addOption(CENTER_CSV, true, "Path towards a CSV containing center (hospital) data.");
        options.addOption(SIGNATURE, true, "Path towards a image file containing the signature to be appended at the end of the report.");
        options.addOption(RVA_LOGO, true, "Path towards a image file containing the logo.");
        options.addOption(COMMENTS, true, "Additional comments to be added to the report, if any.");
        options.addOption(FUSION_CSV, true, "Path towards a CSV of fusions of the sample");
        options.addOption(DISRUPTION_CSV, true, "Path towards a CSV of disruptions of the sample");
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
