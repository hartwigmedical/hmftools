package com.hartwig.hmftools.healthchecker;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.healthchecker.result.QCValue;
import com.hartwig.hmftools.healthchecker.runners.FlagstatChecker;
import com.hartwig.hmftools.healthchecker.runners.HealthChecker;
import com.hartwig.hmftools.healthchecker.runners.MetricsChecker;
import com.hartwig.hmftools.healthchecker.runners.PurpleChecker;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class HealthChecksApplication {

    private static final Logger LOGGER = LogManager.getLogger(HealthChecksApplication.class);

    private static final String REF_SAMPLE = "reference";
    private static final String TUMOR_SAMPLE = "tumor";
    private static final String REF_WGS_METRICS_FILE = "ref_wgs_metrics_file";
    private static final String TUMOR_WGS_METRICS_FILE = "tum_wgs_metrics_file";
    private static final String REF_FLAGSTAT_FILE = "ref_flagstat_file";
    private static final String TUMOR_FLAGSTAT_FILE = "tum_flagstat_file";
    private static final String PURPLE_DIR = "purple_dir";
    private static final String DO_NOT_WRITE_EVALUATION_FILE = "do_not_write_evaluation_file";
    private static final String OUTPUT_DIR = "output_dir";

    @NotNull
    private final String refSample;
    @Nullable
    private final String tumorSample;
    @NotNull
    private final String refWgsMetricsFile;
    @Nullable
    private final String tumorWgsMetricsFile;
    @NotNull
    private final String refFlagstatFile;
    @Nullable
    private final String tumorFlagstatFile;
    @Nullable
    private final String purpleDir;
    @Nullable
    private final String outputDir;

    @VisibleForTesting
    HealthChecksApplication(@NotNull String refSample, @Nullable String tumorSample, @NotNull String refWgsMetricsFile,
            @Nullable String tumorWgsMetricsFile, @NotNull String refFlagstatFile, @Nullable String tumorFlagstatFile,
            @Nullable String purpleDir, @Nullable String outputDir) {
        this.refSample = refSample;
        this.tumorSample = tumorSample;
        this.refWgsMetricsFile = refWgsMetricsFile;
        this.tumorWgsMetricsFile = tumorWgsMetricsFile;
        this.refFlagstatFile = refFlagstatFile;
        this.tumorFlagstatFile = tumorFlagstatFile;
        this.purpleDir = purpleDir;
        this.outputDir = outputDir;
    }

    public static void main(String... args) throws ParseException, IOException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        boolean writeEvaluationFile = !cmd.hasOption(DO_NOT_WRITE_EVALUATION_FILE);
        String outputDir = cmd.hasOption(OUTPUT_DIR) ? cmd.getOptionValue(OUTPUT_DIR) : null;
        String refSample = cmd.getOptionValue(REF_SAMPLE);
        String refFlagstat = cmd.getOptionValue(REF_FLAGSTAT_FILE);
        String refWgsMetricsFile = cmd.getOptionValue(REF_WGS_METRICS_FILE);

        if (refSample == null || refFlagstat == null || refWgsMetricsFile == null || (writeEvaluationFile && outputDir == null)) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Health-Checker", options);
            System.exit(1);
        }

        String tumorSample = cmd.hasOption(TUMOR_SAMPLE) ? cmd.getOptionValue(TUMOR_SAMPLE) : null;
        String tumorWgsMetricsFile = cmd.hasOption(TUMOR_WGS_METRICS_FILE) ? cmd.getOptionValue(TUMOR_WGS_METRICS_FILE) : null;
        String tumorFlagstat = cmd.hasOption(TUMOR_FLAGSTAT_FILE) ? cmd.getOptionValue(TUMOR_FLAGSTAT_FILE) : null;
        String purpleDir = cmd.hasOption(PURPLE_DIR) ? cmd.getOptionValue(PURPLE_DIR) : null;

        new HealthChecksApplication(refSample,
                tumorSample,
                refWgsMetricsFile,
                tumorWgsMetricsFile,
                refFlagstat,
                tumorFlagstat,
                purpleDir,
                outputDir).run(writeEvaluationFile);
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        options.addOption(REF_SAMPLE, true, "The name of the reference sample");
        options.addOption(TUMOR_SAMPLE, true, "The name of the tumor sample");
        options.addOption(PURPLE_DIR, true, "The directory holding the purple output");
        options.addOption(REF_WGS_METRICS_FILE, true, "The path to the wgs metrics file of reference sample");
        options.addOption(TUMOR_WGS_METRICS_FILE, true, "The path to the wgs metrics file of tumor sample");
        options.addOption(REF_FLAGSTAT_FILE, true, "The path to the flagstat file of reference sample");
        options.addOption(TUMOR_FLAGSTAT_FILE, true, "The path to the flagstat file of tumor sample");
        options.addOption(DO_NOT_WRITE_EVALUATION_FILE, false, "Do not write final success or failure file");
        options.addOption(OUTPUT_DIR, true, "The directory where health checker will write output to");
        return options;
    }

    @VisibleForTesting
    void run(boolean writeEvaluationFile) throws IOException {
        List<HealthChecker> checkers;
        if (tumorSample == null || purpleDir == null) {
            LOGGER.info("Running in SingleSample mode");
            checkers = Lists.newArrayList(new MetricsChecker(refWgsMetricsFile, null), new FlagstatChecker(refFlagstatFile, null));
        } else {
            LOGGER.info("Running in Somatic mode");
            checkers = Lists.newArrayList(new MetricsChecker(refWgsMetricsFile, tumorWgsMetricsFile),
                    new FlagstatChecker(refFlagstatFile, tumorFlagstatFile),
                    new PurpleChecker(tumorSample, purpleDir));
        }

        List<QCValue> qcValues = Lists.newArrayList();
        for (HealthChecker checker : checkers) {
            qcValues.addAll(checker.run());
        }

        for (QCValue qcValue : qcValues) {
            LOGGER.info("QC Metric '{}' has value '{}'", qcValue.type(), qcValue.value());
        }

        if (HealthCheckEvaluation.isPass(qcValues)) {
            if (writeEvaluationFile) {
                String evaluationFile = fileOutputBasePath() + ".HealthCheckSucceeded";
                new FileOutputStream(evaluationFile).close();
                LOGGER.info("Evaluation file written to {}", evaluationFile);
            }
            LOGGER.info("Health check evaluation succeeded.");
        } else {
            if (writeEvaluationFile) {
                String evaluationFile = fileOutputBasePath() + ".HealthCheckFailed";
                new FileOutputStream(evaluationFile).close();
                LOGGER.info("Evaluation file written to {}", evaluationFile);
            }
            LOGGER.info("Health check evaluation failed!");
        }
    }

    @NotNull
    private String fileOutputBasePath() {
        assert outputDir != null;

        String sample = tumorSample != null ? tumorSample : refSample;
        return outputDir + File.separator + sample;
    }
}
