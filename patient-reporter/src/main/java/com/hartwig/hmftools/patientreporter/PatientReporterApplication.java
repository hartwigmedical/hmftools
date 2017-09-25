package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.sql.SQLException;

import javax.xml.stream.XMLStreamException;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.formstatus.ImmutableFormStatusModel;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.gene.GeneModel;
import com.hartwig.hmftools.common.lims.LimsJsonModel;
import com.hartwig.hmftools.patientreporter.algo.NotSequenceableReason;
import com.hartwig.hmftools.patientreporter.algo.NotSequenceableReporter;
import com.hartwig.hmftools.patientreporter.algo.NotSequenceableStudy;
import com.hartwig.hmftools.patientreporter.algo.PatientReporter;
import com.hartwig.hmftools.patientreporter.data.COSMICGeneFusionModel;
import com.hartwig.hmftools.patientreporter.report.PDFWriter;
import com.hartwig.hmftools.patientreporter.report.ReportWriter;
import com.hartwig.hmftools.patientreporter.variants.StructuralVariantAnalyzer;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalyzer;
import com.hartwig.hmftools.svannotation.MySQLAnnotator;
import com.hartwig.hmftools.svannotation.NullAnnotator;
import com.hartwig.hmftools.svannotation.VariantAnnotator;

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

    public static final String VERSION = "3.14";

    // TODO (KODU): Make version read from POM
    //    public static final String VERSION = PatientReporterApplication.class.getPackage().getImplementationVersion();

    private static final String CPCT_ECRF = "cpct_ecrf";
    private static final String LIMS_JSON = "lims_json";
    private static final String REPORT_DIRECTORY = "report_dir";
    private static final String RUN_DIRECTORY = "run_dir";
    private static final String NOT_SEQUENCEABLE = "not_sequenceable";
    private static final String NOT_SEQUENCEABLE_REASON = "not_sequenceable_reason";
    private static final String NOT_SEQUENCEABLE_SAMPLE = "not_sequenceable_sample";
    private static final String DRUP_GENES_CSV = "drup_genes_csv";
    private static final String COSMIC_CSV = "cosmic_csv";
    private static final String ENSEMBL_DB = "ensembl_db";
    private static final String FUSION_CSV = "fusion_csv";

    private static final String CENTER_CSV = "center_csv";
    private static final String SIGNATURE = "signature";

    public static void main(final String... args)
            throws ParseException, IOException, HartwigException, DRException, XMLStreamException, SQLException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        if (cmd.hasOption(NOT_SEQUENCEABLE) && validInputForNonSequenceableReport(cmd)) {
            final String notSequenceableSample = cmd.getOptionValue(NOT_SEQUENCEABLE_SAMPLE);
            final NotSequenceableStudy study = NotSequenceableStudy.fromSample(notSequenceableSample);
            if (study == null) {
                LOGGER.warn("Could not determine study for sample " + notSequenceableSample);
            } else {
                LOGGER.info("Generating non-sequenceable report for " + notSequenceableSample);

                final NotSequenceableReason notSequenceableReason =
                        NotSequenceableReason.fromIdentifier(cmd.getOptionValue(NOT_SEQUENCEABLE_REASON));

                final NotSequenceableReporter reporter =
                        new NotSequenceableReporter(buildCpctEcrfModel(cmd), buildLimsModel(cmd), buildReportWriter(cmd));

                reporter.run(notSequenceableSample, notSequenceableReason, study);
            }
        } else if (validInputForPatientReporter(cmd)) {
            LOGGER.info("Running patient reporter v" + VERSION);

            final HmfReporterData reporterData = buildReporterData(cmd);
            final PatientReporter reporter = buildReporter(reporterData.geneModel(), reporterData.fusionModel(), cmd);

            final PatientReport report = reporter.run(cmd.getOptionValue(RUN_DIRECTORY));
            buildReportWriter(cmd).writeSequenceReport(report, reporterData);
        } else {
            printUsageAndExit(options);
        }
    }

    private static HmfReporterData buildReporterData(@NotNull final CommandLine cmd) throws IOException, HartwigException {
        return HmfReporterDataLoader.buildFromFiles(cmd.getOptionValue(DRUP_GENES_CSV), cmd.getOptionValue(COSMIC_CSV),
                cmd.getOptionValue(CENTER_CSV), cmd.getOptionValue(SIGNATURE), cmd.getOptionValue(FUSION_CSV));
    }

    @NotNull
    private static LimsJsonModel buildLimsModel(@NotNull final CommandLine cmd) throws IOException, EmptyFileException {
        return LimsJsonModel.readModelFromFile(cmd.getOptionValue(LIMS_JSON));
    }

    @NotNull
    private static PatientReporter buildReporter(@NotNull final GeneModel geneModel, @NotNull COSMICGeneFusionModel cosmic,
            @NotNull final CommandLine cmd) throws IOException, EmptyFileException, XMLStreamException, SQLException {
        final VariantAnalyzer variantAnalyzer = VariantAnalyzer.fromSlicingRegions(geneModel);

        final VariantAnnotator annotator;
        if (cmd.hasOption(ENSEMBL_DB)) {
            final String url = "jdbc:" + cmd.getOptionValue(ENSEMBL_DB);
            LOGGER.info("connecting to: {}", url);
            annotator = MySQLAnnotator.make(url);
        } else {
            annotator = NullAnnotator.make();
        }
        final StructuralVariantAnalyzer svAnalyzer = new StructuralVariantAnalyzer(annotator, geneModel.hmfRegions(), cosmic);

        return new PatientReporter(buildCpctEcrfModel(cmd), buildLimsModel(cmd), variantAnalyzer, svAnalyzer);
    }

    @NotNull
    private static CpctEcrfModel buildCpctEcrfModel(@NotNull final CommandLine cmd) throws FileNotFoundException, XMLStreamException {
        LOGGER.info(" Loading ECRF database...");
        final CpctEcrfModel cpctEcrfModel =
                CpctEcrfModel.loadFromXML(cmd.getOptionValue(CPCT_ECRF), new ImmutableFormStatusModel(Maps.newHashMap()));
        LOGGER.info("  Loaded data for " + cpctEcrfModel.patientCount() + " patients.");
        return cpctEcrfModel;
    }

    @NotNull
    private static ReportWriter buildReportWriter(@NotNull final CommandLine cmd) {
        return new PDFWriter(cmd.getOptionValue(REPORT_DIRECTORY));
    }

    private static void printUsageAndExit(@NotNull final Options options) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Patient-Reporter", options);
        System.exit(1);
    }

    private static boolean validInputForPatientReporter(@NotNull final CommandLine cmd) {
        if (validInputForEcrfAndTumorPercentages(cmd) && validInputForReportWriter(cmd)) {
            final String drupGenesCsv = cmd.getOptionValue(DRUP_GENES_CSV);
            final String cosmicCsv = cmd.getOptionValue(COSMIC_CSV);
            final String centerCsv = cmd.getOptionValue(CENTER_CSV);
            final String runDirectory = cmd.getOptionValue(RUN_DIRECTORY);
            final String signaturePath = cmd.getOptionValue(SIGNATURE);

            if (drupGenesCsv == null || !exists(drupGenesCsv)) {
                LOGGER.warn(DRUP_GENES_CSV + " has to be an existing file: " + drupGenesCsv);
            } else if (cosmicCsv == null || !exists(cosmicCsv)) {
                LOGGER.warn(COSMIC_CSV + " has to be an existing file: " + cosmicCsv);
            } else if (centerCsv == null || !exists(centerCsv)) {
                LOGGER.warn(CENTER_CSV + " has to be an existing file: " + centerCsv);
            } else if (signaturePath == null || !exists(signaturePath)) {
                LOGGER.warn(SIGNATURE + " has to be an existing file: " + signaturePath);
            } else if (runDirectory == null || !exists(runDirectory) && !isDirectory(runDirectory)) {
                LOGGER.warn(RUN_DIRECTORY + " has to be an existing directory: " + runDirectory);
            } else {
                return true;
            }
        }

        return false;
    }

    private static boolean validInputForNonSequenceableReport(@NotNull final CommandLine cmd) {
        if (validInputForEcrfAndTumorPercentages(cmd) && validInputForReportWriter(cmd)) {
            final NotSequenceableReason notSequenceableReason =
                    NotSequenceableReason.fromIdentifier(cmd.getOptionValue(NOT_SEQUENCEABLE_REASON));
            final String notSequenceableSample = cmd.getOptionValue(NOT_SEQUENCEABLE_SAMPLE);

            if (notSequenceableReason == NotSequenceableReason.OTHER) {
                LOGGER.warn(NOT_SEQUENCEABLE_REASON + " has to be either low_tumor_percentage or low_dna_yield.");
            } else if (notSequenceableSample == null) {
                LOGGER.warn(NOT_SEQUENCEABLE_SAMPLE + " has to be provided.");
            } else {
                return true;
            }
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

    private static boolean validInputForEcrfAndTumorPercentages(@NotNull final CommandLine cmd) {
        final String cpctEcrf = cmd.getOptionValue(CPCT_ECRF);
        final String limsJson = cmd.getOptionValue(LIMS_JSON);

        if (cpctEcrf == null || !exists(cpctEcrf)) {
            LOGGER.warn(CPCT_ECRF + " has to be an existing file: " + cpctEcrf);
        } else if (limsJson == null || !exists(limsJson)) {
            LOGGER.warn(LIMS_JSON + " has to be an existing file: " + limsJson);
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
        options.addOption(CPCT_ECRF, true, "Complete path towards the cpct ecrf xml database.");
        options.addOption(LIMS_JSON, true, "Complete path towards a JSON containing the LIMS data dump.");
        options.addOption(REPORT_DIRECTORY, true, "Complete path to where the PDF reports have to be saved.");
        options.addOption(RUN_DIRECTORY, true, "Complete path towards a single run dir where patient reporter will run on.");
        options.addOption(NOT_SEQUENCEABLE, false, "If set, generates a non-sequenceable report.");
        options.addOption(NOT_SEQUENCEABLE_REASON, true, "Either 'low_tumor_percentage' or 'low_dna_yield'");
        options.addOption(NOT_SEQUENCEABLE_SAMPLE, true, "In case of non-sequenceable reports, the name of the sample used.");
        options.addOption(DRUP_GENES_CSV, true, "Path towards a CSV containing genes that could potentially indicate inclusion in DRUP.");
        options.addOption(COSMIC_CSV, true, "Path towards a CSV containing COSMIC census data.");
        options.addOption(ENSEMBL_DB, true, "Annotate structural variants using this Ensembl DB URI");
        options.addOption(CENTER_CSV, true, "Path towards a CSV containing center data.");
        options.addOption(SIGNATURE, true, "Path towards a image file containing the signature to be appended at the end of the report.");
        options.addOption(FUSION_CSV, true, "Path towards a CSV containing white-listed gene fusions.");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
