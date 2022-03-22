package com.hartwig.hmftools.healthchecker;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.healthchecker.result.QCValue;
import com.hartwig.hmftools.healthchecker.runners.FlagstatChecker;
import com.hartwig.hmftools.healthchecker.runners.HealthCheckSampleConfiguration;
import com.hartwig.hmftools.healthchecker.runners.HealthChecker;
import com.hartwig.hmftools.healthchecker.runners.ReferenceFlagstatChecker;
import com.hartwig.hmftools.healthchecker.runners.ReferenceMetricsChecker;
import com.hartwig.hmftools.healthchecker.runners.TumorFlagstatChecker;
import com.hartwig.hmftools.healthchecker.runners.TumorMetricsChecker;
import com.hartwig.hmftools.healthchecker.runners.PurpleChecker;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class HealthChecksApplication {

    private static final Logger LOGGER = LogManager.getLogger(HealthChecksApplication.class);

    private static final String APPLICATION = "Health-Checker";
    private static final String VERSION = HealthChecksApplication.class.getPackage().getImplementationVersion();

    private static final String REF_SAMPLE = "reference";
    private static final String TUMOR_SAMPLE = "tumor";
    private static final String REF_WGS_METRICS_FILE = "ref_wgs_metrics_file";
    private static final String TUMOR_WGS_METRICS_FILE = "tum_wgs_metrics_file";
    private static final String REF_FLAGSTAT_FILE = "ref_flagstat_file";
    private static final String TUMOR_FLAGSTAT_FILE = "tum_flagstat_file";
    private static final String PURPLE_DIR = "purple_dir";
    private static final String DO_NOT_WRITE_EVALUATION_FILE = "do_not_write_evaluation_file";
    private static final String OUTPUT_DIR = "output_dir";

    @Nullable
    private final HealthCheckSampleConfiguration refSampleConfig;
    @Nullable
    private final HealthCheckSampleConfiguration tumorSampleConfig;
    @Nullable
    private final String purpleDir;
    @Nullable
    private final String outputDir;

    @VisibleForTesting
    public HealthChecksApplication(@Nullable final HealthCheckSampleConfiguration refSampleConfig,
            @Nullable final HealthCheckSampleConfiguration tumorSampleConfig, @Nullable final String purpleDir,
            @Nullable final String outputDir) {
        this.refSampleConfig = refSampleConfig;
        this.tumorSampleConfig = tumorSampleConfig;
        this.purpleDir = purpleDir;
        this.outputDir = outputDir;
    }

    public static void main(String... args) throws ParseException, IOException {
        LOGGER.info("Running {} v{}", APPLICATION, VERSION);
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        boolean writeEvaluationFile = !cmd.hasOption(DO_NOT_WRITE_EVALUATION_FILE);
        String outputDir = cmd.hasOption(OUTPUT_DIR) ? cmd.getOptionValue(OUTPUT_DIR) : null;
        if (writeEvaluationFile && outputDir == null) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp(APPLICATION, options);
            System.exit(1);
        }

        String refSample = cmd.getOptionValue(REF_SAMPLE, null);
        String refFlagstat = cmd.getOptionValue(REF_FLAGSTAT_FILE, null);
        String refWgsMetricsFile = cmd.getOptionValue(REF_WGS_METRICS_FILE, null);
        String tumorSample = cmd.getOptionValue(TUMOR_SAMPLE, null);
        String tumorWgsMetricsFile = cmd.getOptionValue(TUMOR_WGS_METRICS_FILE, null);
        String tumorFlagstat = cmd.getOptionValue(TUMOR_FLAGSTAT_FILE, null);

        String purpleDir = cmd.getOptionValue(PURPLE_DIR, null);
        new HealthChecksApplication(HealthCheckSampleConfiguration.of(tumorSample, tumorWgsMetricsFile, tumorFlagstat),
                HealthCheckSampleConfiguration.of(refSample, refWgsMetricsFile, refWgsMetricsFile),
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
        Optional<HealthCheckSampleConfiguration> maybeRefSampleConfiguration = Optional.ofNullable(refSampleConfig);
        Optional<HealthCheckSampleConfiguration> maybeTumorSampleConfiguration = Optional.ofNullable(tumorSampleConfig);

        if (maybeRefSampleConfiguration.isPresent() && maybeTumorSampleConfiguration.isEmpty()) {
            LOGGER.info("Running in germline only mode");
            checkers = List.of(new ReferenceMetricsChecker(maybeRefSampleConfiguration.get().wgsMetricsFile()),
                    new ReferenceFlagstatChecker(maybeRefSampleConfiguration.get().flagstatFile()));
        } else if (maybeRefSampleConfiguration.isEmpty() && maybeTumorSampleConfiguration.isPresent()) {
            LOGGER.info("Running in tumor only mode");
            checkers = List.of(new TumorMetricsChecker(maybeTumorSampleConfiguration.get().wgsMetricsFile()),
                    new TumorFlagstatChecker(maybeTumorSampleConfiguration.get().flagstatFile()));
        } else if (maybeTumorSampleConfiguration.isPresent() && purpleDir != null) {
            LOGGER.info("Running in somatic mode");
            checkers = Lists.newArrayList(new TumorMetricsChecker(maybeTumorSampleConfiguration.get().wgsMetricsFile()),
                    new ReferenceMetricsChecker(maybeRefSampleConfiguration.get().wgsMetricsFile()),
                    new TumorFlagstatChecker(maybeTumorSampleConfiguration.get().flagstatFile()),
                    new ReferenceFlagstatChecker(maybeRefSampleConfiguration.get().flagstatFile()),
                    new PurpleChecker(maybeTumorSampleConfiguration.get().sampleName(), purpleDir));
        } else {
            throw new IllegalArgumentException(String.format("Illegal combination of arguments: [%s, %s, %s]",
                    maybeRefSampleConfiguration,
                    maybeTumorSampleConfiguration,
                    purpleDir));
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
        Optional<String> tumorSample = Optional.ofNullable(tumorSampleConfig).map(HealthCheckSampleConfiguration::sampleName);
        Optional<String> refSample = Optional.ofNullable(refSampleConfig).map(HealthCheckSampleConfiguration::sampleName);
        return Optional.ofNullable(outputDir)
                .map(o -> outputDir + File.separator + tumorSample.orElseGet(refSample::orElseThrow))
                .orElseThrow();
    }
}
