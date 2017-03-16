package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import javax.xml.stream.XMLStreamException;

import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.slicing.Slicer;
import com.hartwig.hmftools.common.slicing.SlicerFactory;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalyzer;
import com.hartwig.hmftools.patientreporter.lims.TumorPercentages;
import com.hartwig.hmftools.patientreporter.report.PDFWriter;
import com.hartwig.hmftools.patientreporter.report.ReportWriter;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalyzer;

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

    private static final String CPCT_SLICING_BED_ARGS_DESC = "Complete path towards the CPCT slicing bed.";
    private static final String CPCT_SLICING_BED = "cpct_slicing_bed";

    private static final String HIGH_CONFIDENCE_BED_ARGS_DESC = "Complete path towards the high confidence bed.";
    private static final String HIGH_CONFIDENCE_BED = "high_confidence_bed";

    private static final String HMF_SLICING_BED_ARGS_DESC = "Complete path towards the HMF slicing bed.";
    private static final String HMF_SLICING_BED = "hmf_slicing_bed";

    private static final String CPCT_ECRF_ARGS_DESC = "Complete path towards the cpct ecrf xml database.";
    private static final String CPCT_ECRF = "cpct_ecrf";

    private static final String TUMOR_PERCENTAGE_CSV_ARGS_DESC = "Complete path towards a CSV containing tumor percentages.";
    private static final String TUMOR_PERCENTAGE_CSV = "tumor_percentage_csv";

    private static final String REPORT_LOGO_ARGS_DESC = "Complete path to the logo used in the PDF report.";
    private static final String REPORT_LOGO = "report_logo";

    private static final String REPORT_DIRECTORY_ARGS_DESC = "Complete path to where the PDF reports have to be saved.";
    private static final String REPORT_DIRECTORY = "report_dir";

    private static final String TMP_DIRECTORY_ARGS_DESC = "Complete path where, if provided, temporary output files will be written to.";
    private static final String TMP_DIRECTORY = "tmp_dir";

    private static final String RUN_DIRECTORY_ARGS_DESC = "Complete path towards a single rundir where patient reporter will run on.";
    private static final String RUN_DIRECTORY = "run_dir";

    private static final String BATCH_MODE_ARGS_DESC = "If set, runs in batch mode and generates a CSV with statistics.";
    private static final String BATCH_MODE = "batch_mode";

    private static final String BATCH_DIRECTORY_ARGS_DESC = "The directory that will be iterated over in batch-mode.";
    private static final String BATCH_DIRECTORY = "batch_directory";

    private static final String BATCH_OUTPUT_ARGS_DESC = "The file which will contain the results of the batch mode output.";
    private static final String BATCH_OUTPUT = "batch_output";

    private static final String NOT_SEQUENCEABLE_ARGS_DESC = "If set, generates a non-sequenceable report.";
    private static final String NOT_SEQUENCEABLE = "not_sequenceable";

    private static final String NOT_SEQUENCEABLE_REASON_ARGS_DESC = "Either 'low_tumor_percentage' or 'low_dna_yield'";
    private static final String NOT_SEQUENCEABLE_REASON = "not_sequenceable_reason";

    private static final String NOT_SEQUENCEABLE_SAMPLE_ARGS_DESC = "In case of non-sequenceable reports, the name of the sample used.";
    private static final String NOT_SEQUENCEABLE_SAMPLE = "not_sequenceable_sample";

    public static void main(final String... args)
            throws ParseException, IOException, HartwigException, DRException, XMLStreamException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        if (cmd.hasOption(NOT_SEQUENCEABLE) && validInputForNonSequenceableReport(cmd)) {
            final String notSequenceableSample = cmd.getOptionValue(NOT_SEQUENCEABLE_SAMPLE);
            LOGGER.info("Generating non-sequenceable report for " + notSequenceableSample);

            final TumorPercentages tumorPercentages = buildTumorPercentages(cmd);
            final ReportWriter reportWriter = buildReportWriter(cmd);

            final NotSequenceableReason notSequenceableReason = NotSequenceableReason.fromIdentifier(
                    cmd.getOptionValue(NOT_SEQUENCEABLE_REASON));

            final NotSequenceableReporter reporter = new NotSequenceableReporter(buildCpctEcrfModel(cmd), reportWriter,
                    tumorPercentages);

            reporter.run(notSequenceableSample, notSequenceableReason);
        } else if (validInputForPatientReporter(cmd)) {
            LOGGER.info("Running patient reporter");
            final Slicer hmfSlicingRegion = buildHmfSlicingRegion(cmd);
            final SinglePatientReporter reporter = buildReporter(hmfSlicingRegion, cmd);

            if (cmd.hasOption(BATCH_MODE) && validInputForBatchMode(cmd)) {
                LOGGER.info("Switching to running patient reporter in batch-mode.");
                final BatchReportAnalyser analyser = new BatchReportAnalyser(reporter);
                final List<String> batchAnalysis = analyser.run(cmd.getOptionValue(BATCH_DIRECTORY));
                Files.write(new File(cmd.getOptionValue(BATCH_OUTPUT)).toPath(), batchAnalysis);
            } else if (validInputForSinglePatientReport(cmd)) {
                final PatientReport report = reporter.run(cmd.getOptionValue(RUN_DIRECTORY));
                buildReportWriter(cmd).writeSequenceReport(report, hmfSlicingRegion);
            } else {
                printUsageAndExit(options);
            }
        } else {
            printUsageAndExit(options);
        }
    }

    @NotNull
    private static Slicer buildHmfSlicingRegion(@NotNull final CommandLine cmd)
            throws IOException, EmptyFileException {
        return SlicerFactory.fromBedFile(cmd.getOptionValue(HMF_SLICING_BED));
    }

    @NotNull
    private static TumorPercentages buildTumorPercentages(@NotNull final CommandLine cmd)
            throws IOException, EmptyFileException {
        return TumorPercentages.loadFromCsv(cmd.getOptionValue(TUMOR_PERCENTAGE_CSV));
    }

    @NotNull
    private static SinglePatientReporter buildReporter(@NotNull final Slicer hmfSlicingRegion,
            @NotNull final CommandLine cmd) throws IOException, EmptyFileException, XMLStreamException {
        final VariantAnalyzer variantAnalyzer = VariantAnalyzer.fromSlicingRegions(hmfSlicingRegion,
                SlicerFactory.fromBedFile(cmd.getOptionValue(HIGH_CONFIDENCE_BED)),
                SlicerFactory.fromBedFile(cmd.getOptionValue(CPCT_SLICING_BED)));
        final CopyNumberAnalyzer copyNumberAnalyzer = CopyNumberAnalyzer.fromHmfSlicingRegion(hmfSlicingRegion);

        return new SinglePatientReporter(buildCpctEcrfModel(cmd), variantAnalyzer, copyNumberAnalyzer,
                buildTumorPercentages(cmd), cmd.getOptionValue(TMP_DIRECTORY));
    }

    @NotNull
    private static CpctEcrfModel buildCpctEcrfModel(@NotNull final CommandLine cmd)
            throws FileNotFoundException, XMLStreamException {
        LOGGER.info(" Loading ECRF database...");
        final CpctEcrfModel cpctEcrfModel = CpctEcrfModel.loadFromXML(cmd.getOptionValue(CPCT_ECRF));
        LOGGER.info("  Loaded data for " + cpctEcrfModel.patientCount() + " patients.");
        return cpctEcrfModel;
    }

    @NotNull
    private static ReportWriter buildReportWriter(@NotNull final CommandLine cmd) {
        return new PDFWriter(cmd.getOptionValue(REPORT_DIRECTORY), cmd.getOptionValue(REPORT_LOGO));
    }

    private static void printUsageAndExit(@NotNull final Options options) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Patient-Reporter", options);
        System.exit(1);
    }

    private static boolean validInputForSinglePatientReport(@NotNull final CommandLine cmd) {
        if (validInputForPatientReporter(cmd) && validInputForReportWriter(cmd)) {
            final String runDirectory = cmd.getOptionValue(RUN_DIRECTORY);

            if (runDirectory == null || !exists(runDirectory) && !isDirectory(runDirectory)) {
                LOGGER.warn(RUN_DIRECTORY + " has to be an existing directory: " + runDirectory);
            } else {
                return true;
            }
        }

        return false;
    }

    private static boolean validInputForBatchMode(@NotNull final CommandLine cmd) {
        if (validInputForPatientReporter(cmd)) {
            final String batchDirectory = cmd.getOptionValue(BATCH_DIRECTORY);
            final String batchOutputFile = cmd.getOptionValue(BATCH_OUTPUT);

            if (batchDirectory == null || !exists(batchDirectory) && !isDirectory(batchDirectory)) {
                LOGGER.warn(BATCH_DIRECTORY + " has to be an existing directory: " + batchDirectory);
            } else if (batchOutputFile == null) {
                LOGGER.warn(BATCH_OUTPUT + " has to be provided.");
            } else {
                return true;
            }
        }

        return false;
    }

    private static boolean validInputForPatientReporter(@NotNull final CommandLine cmd) {
        if (validInputForEcrfAndTumorPercentages(cmd)) {
            final String cpctSlicingBed = cmd.getOptionValue(CPCT_SLICING_BED);
            final String highConfidenceBed = cmd.getOptionValue(HIGH_CONFIDENCE_BED);
            final String hmfSlicingBed = cmd.getOptionValue(HMF_SLICING_BED);
            final String tmpDirectory = cmd.getOptionValue(TMP_DIRECTORY);

            if (cpctSlicingBed == null || !exists(cpctSlicingBed)) {
                LOGGER.warn(CPCT_SLICING_BED + " has to be an existing file: " + cpctSlicingBed);
            } else if (highConfidenceBed == null || !exists(highConfidenceBed)) {
                LOGGER.warn(HIGH_CONFIDENCE_BED + " has to be an existing file: " + highConfidenceBed);
            } else if (hmfSlicingBed == null || !exists(hmfSlicingBed)) {
                LOGGER.warn(HMF_SLICING_BED + " has to be an existing file: " + hmfSlicingBed);
            } else if (tmpDirectory != null && (!exists(tmpDirectory) || !isDirectory(tmpDirectory))) {
                LOGGER.warn(TMP_DIRECTORY + " has to be an existing directory: " + highConfidenceBed);
            } else {
                return true;
            }
        }

        return false;
    }

    private static boolean validInputForNonSequenceableReport(@NotNull final CommandLine cmd) {
        if (validInputForEcrfAndTumorPercentages(cmd) && validInputForReportWriter(cmd)) {
            final NotSequenceableReason notSequenceableReason = NotSequenceableReason.fromIdentifier(
                    cmd.getOptionValue(NOT_SEQUENCEABLE_REASON));
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
        final String reportLogo = cmd.getOptionValue(REPORT_LOGO);
        final String reportDirectory = cmd.getOptionValue(REPORT_DIRECTORY);

        if (reportLogo == null || !exists(reportLogo)) {
            LOGGER.warn(REPORT_LOGO + " has to be an existing file: " + reportLogo);
        } else if (reportDirectory == null || !exists(reportDirectory) || !isDirectory(reportDirectory)) {
            LOGGER.warn(REPORT_DIRECTORY + " has to be an existing directory: " + reportDirectory);
        } else {
            return true;
        }

        return false;
    }

    private static boolean validInputForEcrfAndTumorPercentages(@NotNull final CommandLine cmd) {
        final String cpctEcrf = cmd.getOptionValue(CPCT_ECRF);
        final String tumorPercentageCsv = cmd.getOptionValue(TUMOR_PERCENTAGE_CSV);

        if (cpctEcrf == null || !exists(cpctEcrf)) {
            LOGGER.warn(CPCT_ECRF + " has to be an existing file: " + cpctEcrf);
        } else if (tumorPercentageCsv == null || !exists(tumorPercentageCsv)) {
            LOGGER.warn(TUMOR_PERCENTAGE_CSV + " has to be an existing file: " + tumorPercentageCsv);
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

        options.addOption(CPCT_SLICING_BED, true, CPCT_SLICING_BED_ARGS_DESC);
        options.addOption(HIGH_CONFIDENCE_BED, true, HIGH_CONFIDENCE_BED_ARGS_DESC);
        options.addOption(HMF_SLICING_BED, true, HMF_SLICING_BED_ARGS_DESC);
        options.addOption(CPCT_ECRF, true, CPCT_ECRF_ARGS_DESC);
        options.addOption(TUMOR_PERCENTAGE_CSV, true, TUMOR_PERCENTAGE_CSV_ARGS_DESC);
        options.addOption(REPORT_LOGO, true, REPORT_LOGO_ARGS_DESC);
        options.addOption(REPORT_DIRECTORY, true, REPORT_DIRECTORY_ARGS_DESC);
        options.addOption(TMP_DIRECTORY, true, TMP_DIRECTORY_ARGS_DESC);
        options.addOption(RUN_DIRECTORY, true, RUN_DIRECTORY_ARGS_DESC);

        options.addOption(BATCH_MODE, false, BATCH_MODE_ARGS_DESC);
        options.addOption(BATCH_DIRECTORY, true, BATCH_DIRECTORY_ARGS_DESC);
        options.addOption(BATCH_OUTPUT, true, BATCH_OUTPUT_ARGS_DESC);

        options.addOption(NOT_SEQUENCEABLE, false, NOT_SEQUENCEABLE_ARGS_DESC);
        options.addOption(NOT_SEQUENCEABLE_REASON, true, NOT_SEQUENCEABLE_REASON_ARGS_DESC);
        options.addOption(NOT_SEQUENCEABLE_SAMPLE, true, NOT_SEQUENCEABLE_SAMPLE_ARGS_DESC);

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
