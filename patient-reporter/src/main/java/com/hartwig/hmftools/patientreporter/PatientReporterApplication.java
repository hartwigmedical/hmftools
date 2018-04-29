package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.common.center.Center;
import com.hartwig.hmftools.common.center.CenterModel;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
import com.hartwig.hmftools.patientreporter.algo.ImmutableNotSequenceableReporter;
import com.hartwig.hmftools.patientreporter.algo.ImmutablePatientReporter;
import com.hartwig.hmftools.patientreporter.algo.NotSequenceableReason;
import com.hartwig.hmftools.patientreporter.algo.NotSequenceableReporter;
import com.hartwig.hmftools.patientreporter.algo.NotSequenceableStudy;
import com.hartwig.hmftools.patientreporter.algo.PatientReporter;
import com.hartwig.hmftools.patientreporter.civic.CivicAnalyzer;
import com.hartwig.hmftools.patientreporter.report.PDFWriter;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalyzer;
import com.hartwig.hmftools.svannotation.MySQLAnnotator;
import com.hartwig.hmftools.svannotation.NullAnnotator;
import com.hartwig.hmftools.svannotation.VariantAnnotator;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalyzer;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.exception.DRException;

public class PatientReporterApplication {

    private static final Logger LOGGER = LogManager.getLogger(PatientReporterApplication.class);

    public static final String VERSION = PatientReporterApplication.class.getPackage().getImplementationVersion();

    // KODU: Uncomment this line when generating an example report using PDFWriterTest
    //    public static final String VERSION = "4.13";

    private static final String TUMOR_LOCATION_CSV = "tumor_location_csv";
    private static final String LIMS_JSON = "lims_json";
    private static final String REPORT_DIRECTORY = "report_dir";
    private static final String RUN_DIRECTORY = "run_dir";
    private static final String NOT_SEQUENCEABLE = "not_sequenceable";
    private static final String NOT_SEQUENCEABLE_REASON = "not_sequenceable_reason";
    private static final String NOT_SEQUENCEABLE_SAMPLE = "not_sequenceable_sample";

    private static final String COSMIC_GENE_CSV = "cosmic_gene_csv";
    private static final String FUSION_PAIRS_CSV = "fusion_pairs_csv";
    private static final String PROMISCUOUS_FIVE_CSV = "promiscuous_five_csv";
    private static final String PROMISCUOUS_THREE_CSV = "promiscuous_three_csv";
    private static final String DRUP_GENES_CSV = "drup_genes_csv";
    private static final String ENSEMBL_DB = "ensembl_db";
    private static final String FASTA_FILE_LOCATION = "fasta_file_location";
    private static final String CENTER_CSV = "center_csv";

    private static final String SIGNATURE = "signature";
    private static final String COMMENTS = "comments";

    public static void main(final String... args) throws ParseException, IOException, DRException, SQLException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        if (!validInputForReportWriter(cmd) || !validInputForBaseReporterData(cmd)) {
            printUsageAndExit(options);
        }
        LOGGER.info("Running patient reporter v" + VERSION);
        final PDFWriter pdfWriter = new PDFWriter(cmd.getOptionValue(REPORT_DIRECTORY));

        if (cmd.hasOption(NOT_SEQUENCEABLE) && validInputForNonSequenceableReport(cmd)) {
            final String sample = cmd.getOptionValue(NOT_SEQUENCEABLE_SAMPLE);
            LOGGER.info("Generating non-sequenceable report for {}", sample);
            final NotSequenceableReason reason = NotSequenceableReason.fromIdentifier(cmd.getOptionValue(NOT_SEQUENCEABLE_REASON));
            final NotSequenceableReporter reporter = ImmutableNotSequenceableReporter.of(buildBaseReporterData(cmd));

            final NotSequencedPatientReport report = reporter.run(sample, reason, cmd.getOptionValue(COMMENTS));
            pdfWriter.writeNonSequenceableReport(report);
        } else if (validInputForPatientReporter(cmd)) {
            LOGGER.info("Generating sequenceable report...");
            final HmfReporterData reporterData = buildReporterData(cmd);
            final PatientReporter reporter = buildReporter(cmd, reporterData);

            final SequencedPatientReport report = reporter.run(cmd.getOptionValue(RUN_DIRECTORY), cmd.getOptionValue(COMMENTS));
            pdfWriter.writeSequenceReport(report, reporterData);
        } else {
            printUsageAndExit(options);
        }
    }

    @NotNull
    private static BaseReporterData buildBaseReporterData(@NotNull final CommandLine cmd) throws IOException {
        LOGGER.info(" Loading ECRF CSV dump...");
        final List<PatientTumorLocation> patientTumorLocations = PatientTumorLocation.readRecords(cmd.getOptionValue(TUMOR_LOCATION_CSV));
        LOGGER.info("  Loaded data for {} patients.", patientTumorLocations.size());
        LOGGER.info(" Loading LIMS database...");
        final Lims lims = LimsFactory.fromLimsJson(cmd.getOptionValue(LIMS_JSON));
        LOGGER.info("  Loaded data for {} samples.", lims.sampleCount());
        final CenterModel centerModel = Center.readFromCSV(cmd.getOptionValue(CENTER_CSV));
        return ImmutableBaseReporterData.of(patientTumorLocations, lims, centerModel, cmd.getOptionValue(SIGNATURE));
    }

    @NotNull
    private static HmfReporterData buildReporterData(@NotNull final CommandLine cmd) throws IOException {
        return HmfReporterDataLoader.buildFromFiles(cmd.getOptionValue(COSMIC_GENE_CSV),
                cmd.getOptionValue(FUSION_PAIRS_CSV),
                cmd.getOptionValue(PROMISCUOUS_FIVE_CSV),
                cmd.getOptionValue(PROMISCUOUS_THREE_CSV),
                cmd.getOptionValue(DRUP_GENES_CSV),
                cmd.getOptionValue(FASTA_FILE_LOCATION));
    }

    @NotNull
    private static PatientReporter buildReporter(@NotNull final CommandLine cmd, @NotNull final HmfReporterData reporterData)
            throws IOException, SQLException {
        final VariantAnalyzer variantAnalyzer = VariantAnalyzer.of(reporterData.panelGeneModel(), reporterData.microsatelliteAnalyzer());

        final VariantAnnotator annotator;
        if (cmd.hasOption(ENSEMBL_DB)) {
            final String url = "jdbc:" + cmd.getOptionValue(ENSEMBL_DB);
            LOGGER.info("connecting to: {}", url);
            annotator = MySQLAnnotator.make(url);
        } else {
            annotator = NullAnnotator.make();
        }
        final StructuralVariantAnalyzer svAnalyzer =
                new StructuralVariantAnalyzer(annotator, reporterData.panelGeneModel().regions(), reporterData.knownFusionsModel());

        return ImmutablePatientReporter.of(buildBaseReporterData(cmd), reporterData, variantAnalyzer, svAnalyzer, new CivicAnalyzer());
    }

    private static boolean validInputForPatientReporter(@NotNull final CommandLine cmd) {
        final String runDirectory = cmd.getOptionValue(RUN_DIRECTORY);
        final String drupGenesCsv = cmd.getOptionValue(DRUP_GENES_CSV);
        final String cosmicGeneCsv = cmd.getOptionValue(COSMIC_GENE_CSV);
        final String fusionPairsCsv = cmd.getOptionValue(FUSION_PAIRS_CSV);
        final String promiscuousFiveCsv = cmd.getOptionValue(PROMISCUOUS_FIVE_CSV);
        final String promiscuousThreeCsv = cmd.getOptionValue(PROMISCUOUS_THREE_CSV);

        final String fastaFileLocation = cmd.getOptionValue(FASTA_FILE_LOCATION);

        if (runDirectory == null || !exists(runDirectory) && !isDirectory(runDirectory)) {
            LOGGER.warn(RUN_DIRECTORY + " has to be an existing directory: " + runDirectory);
        } else if (drupGenesCsv == null || !exists(drupGenesCsv)) {
            LOGGER.warn(DRUP_GENES_CSV + " has to be an existing file: " + drupGenesCsv);
        } else if (cosmicGeneCsv == null || !exists(cosmicGeneCsv)) {
            LOGGER.warn(COSMIC_GENE_CSV + " has to be an existing file: " + cosmicGeneCsv);
        } else if (fusionPairsCsv == null || !exists(fusionPairsCsv)) {
            LOGGER.warn(FUSION_PAIRS_CSV + " has to be an existing file: " + fusionPairsCsv);
        } else if (promiscuousFiveCsv == null || !exists(promiscuousFiveCsv)) {
            LOGGER.warn(PROMISCUOUS_FIVE_CSV + " has to be an existing file: " + promiscuousFiveCsv);
        } else if (promiscuousThreeCsv == null || !exists(promiscuousThreeCsv)) {
            LOGGER.warn(PROMISCUOUS_THREE_CSV + " has to be an existing file: " + promiscuousThreeCsv);
        } else if (fastaFileLocation == null || !exists(fastaFileLocation)) {
            LOGGER.warn(FASTA_FILE_LOCATION + " has to be an existing file: " + fastaFileLocation);
        } else {
            return true;
        }
        return false;
    }

    private static boolean validInputForNonSequenceableReport(@NotNull final CommandLine cmd) {
        final NotSequenceableReason notSequenceableReason =
                NotSequenceableReason.fromIdentifier(cmd.getOptionValue(NOT_SEQUENCEABLE_REASON));
        final String notSequenceableSample = cmd.getOptionValue(NOT_SEQUENCEABLE_SAMPLE);

        if (notSequenceableReason == NotSequenceableReason.UNDEFINED) {
            LOGGER.warn(NOT_SEQUENCEABLE_REASON + " has to be low_tumor_percentage, low_dna_yield or post_isolation_fail.");
        } else if (notSequenceableSample == null) {
            LOGGER.warn(NOT_SEQUENCEABLE_SAMPLE + " has to be provided.");
        } else if (NotSequenceableStudy.fromSample(notSequenceableSample) == null) {
            LOGGER.warn("Could not determine study for sample " + notSequenceableSample);
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

    private static boolean validInputForBaseReporterData(@NotNull final CommandLine cmd) {
        final String tumorLocationCsv = cmd.getOptionValue(TUMOR_LOCATION_CSV);
        final String limsJson = cmd.getOptionValue(LIMS_JSON);
        final String centerCsv = cmd.getOptionValue(CENTER_CSV);
        final String signaturePath = cmd.getOptionValue(SIGNATURE);

        if (tumorLocationCsv == null || !exists(tumorLocationCsv)) {
            LOGGER.warn(TUMOR_LOCATION_CSV + " has to be an existing file: " + tumorLocationCsv);
        } else if (limsJson == null || !exists(limsJson)) {
            LOGGER.warn(LIMS_JSON + " has to be an existing file: " + limsJson);
        } else if (centerCsv == null || !exists(centerCsv)) {
            LOGGER.warn(CENTER_CSV + " has to be an existing file: " + centerCsv);
        } else if (signaturePath == null || !exists(signaturePath)) {
            LOGGER.warn(SIGNATURE + " has to be an existing file: " + signaturePath);
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
        options.addOption(LIMS_JSON, true, "Complete path towards a JSON containing the LIMS data dump.");
        options.addOption(REPORT_DIRECTORY, true, "Complete path to where the PDF reports have to be saved.");
        options.addOption(RUN_DIRECTORY, true, "Complete path towards a single run dir where patient reporter will run on.");
        options.addOption(NOT_SEQUENCEABLE, false, "If set, generates a non-sequenceable report.");
        options.addOption(NOT_SEQUENCEABLE_REASON, true, "Either 'low_tumor_percentage' or 'low_dna_yield'");
        options.addOption(NOT_SEQUENCEABLE_SAMPLE, true, "In case of non-sequenceable reports, the name of the sample used.");
        options.addOption(COSMIC_GENE_CSV, true, "Path towards a CSV containing COSMIC gene data.");
        options.addOption(FUSION_PAIRS_CSV, true, "Path towards a CSV containing white-listed gene fusion pairs.");
        options.addOption(PROMISCUOUS_FIVE_CSV, true, "Path towards a CSV containing white-listed promiscuous 5' genes.");
        options.addOption(PROMISCUOUS_THREE_CSV, true, "Path towards a CSV containing white-listed promiscuous 3' genes.");
        options.addOption(DRUP_GENES_CSV, true, "Path towards a CSV containing genes that could potentially indicate inclusion in DRUP.");
        options.addOption(ENSEMBL_DB, true, "Annotate structural variants using this Ensembl DB URI");
        options.addOption(CENTER_CSV, true, "Path towards a CSV containing center data.");
        options.addOption(SIGNATURE, true, "Path towards a image file containing the signature to be appended at the end of the report.");
        options.addOption(FASTA_FILE_LOCATION, true, "Path towards the FASTA file containing the ref genome.");
        options.addOption(COMMENTS, true, "Additional comments to be added to the report, if any.");
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
