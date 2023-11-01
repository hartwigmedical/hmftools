package com.hartwig.hmftools.healthchecker;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.healthchecker.result.QCValueType.REF_PROPORTION_DUPLICATE;
import static com.hartwig.hmftools.healthchecker.result.QCValueType.REF_PROPORTION_MAPPED;
import static com.hartwig.hmftools.healthchecker.result.QCValueType.TUM_PROPORTION_DUPLICATE;
import static com.hartwig.hmftools.healthchecker.result.QCValueType.TUM_PROPORTION_MAPPED;

import java.io.FileOutputStream;
import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.healthchecker.result.QCValue;
import com.hartwig.hmftools.healthchecker.runners.HealthCheckSampleConfiguration;
import com.hartwig.hmftools.healthchecker.runners.HealthChecker;
import com.hartwig.hmftools.healthchecker.runners.ReferenceMetricsChecker;
import com.hartwig.hmftools.healthchecker.runners.FlagstatChecker;
import com.hartwig.hmftools.healthchecker.runners.TumorMetricsChecker;
import com.hartwig.hmftools.healthchecker.runners.PurpleChecker;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

public class HealthChecksApplication
{
    public static final Logger HC_LOGGER = LogManager.getLogger(HealthChecksApplication.class);

    private static final String REF_SAMPLE = "reference";
    private static final String TUMOR_SAMPLE = "tumor";
    private static final String REF_WGS_METRICS_FILE = "ref_wgs_metrics_file";
    private static final String TUMOR_WGS_METRICS_FILE = "tum_wgs_metrics_file";
    private static final String REF_FLAGSTAT_FILE = "ref_flagstat_file";
    private static final String TUMOR_FLAGSTAT_FILE = "tum_flagstat_file";
    private static final String DO_NOT_WRITE_EVALUATION_FILE = "do_not_write_evaluation_file";

    @Nullable
    private final HealthCheckSampleConfiguration mRefSampleConfig;
    @Nullable
    private final HealthCheckSampleConfiguration mTumorSampleConfig;
    @Nullable
    private final String mPurpleDir;
    @Nullable
    private final String mOutputDir;

    private final boolean mWriteEvaluationFile;

    public HealthChecksApplication(final ConfigBuilder configBuilder)
    {
        String refSample = configBuilder.getValue(REF_SAMPLE);
        String refFlagstat = configBuilder.getValue(REF_FLAGSTAT_FILE);
        String refWgsMetricsFile = configBuilder.getValue(REF_WGS_METRICS_FILE);
        String tumorSample = configBuilder.getValue(TUMOR_SAMPLE);
        String tumorWgsMetricsFile = configBuilder.getValue(TUMOR_WGS_METRICS_FILE);
        String tumorFlagstat = configBuilder.getValue(TUMOR_FLAGSTAT_FILE);

        mRefSampleConfig = refSample != null ? new HealthCheckSampleConfiguration(refSample, refWgsMetricsFile, refFlagstat) : null;
        mTumorSampleConfig = tumorSample != null ? new HealthCheckSampleConfiguration(tumorSample, tumorWgsMetricsFile, tumorFlagstat) : null;
        mPurpleDir = configBuilder.getValue(PURPLE_DIR_CFG);
        mOutputDir = parseOutputDir(configBuilder);
        mWriteEvaluationFile = !configBuilder.hasFlag(DO_NOT_WRITE_EVALUATION_FILE);
    }

    @VisibleForTesting
    public HealthChecksApplication(
            @Nullable final HealthCheckSampleConfiguration refSampleConfig,
            @Nullable final HealthCheckSampleConfiguration tumorSampleConfig, @Nullable final String purpleDir,
            @Nullable final String outputDir)
    {
        mRefSampleConfig = refSampleConfig;
        mTumorSampleConfig = tumorSampleConfig;
        mPurpleDir = purpleDir;
        mOutputDir = outputDir;
        mWriteEvaluationFile = false;
    }

    @VisibleForTesting
    void run() throws IOException
    {
        List<HealthChecker> checkers;
        Optional<HealthCheckSampleConfiguration> maybeRefSampleConfiguration = Optional.ofNullable(mRefSampleConfig);
        Optional<HealthCheckSampleConfiguration> maybeTumorSampleConfiguration = Optional.ofNullable(mTumorSampleConfig);

        if(maybeRefSampleConfiguration.isPresent() && maybeTumorSampleConfiguration.isEmpty())
        {
            HC_LOGGER.info("Running in germline-only mode");

            checkers = List.of(
                    new ReferenceMetricsChecker(maybeRefSampleConfiguration.get().WgsMetricsFile),
                    new FlagstatChecker(maybeRefSampleConfiguration.get().FlagstatFile, REF_PROPORTION_MAPPED, REF_PROPORTION_DUPLICATE));
        }
        else if(maybeRefSampleConfiguration.isEmpty() && maybeTumorSampleConfiguration.isPresent())
        {
            HC_LOGGER.info("Running in tumor-only mode");

            checkers = List.of(
                    new TumorMetricsChecker(maybeTumorSampleConfiguration.get().WgsMetricsFile),
                    new FlagstatChecker(maybeTumorSampleConfiguration.get().FlagstatFile, TUM_PROPORTION_MAPPED, TUM_PROPORTION_DUPLICATE));
        }
        else if(maybeTumorSampleConfiguration.isPresent() && mPurpleDir != null)
        {
            HC_LOGGER.info("Running in tumor-normal mode");

            checkers = Lists.newArrayList(new TumorMetricsChecker(maybeTumorSampleConfiguration.get().WgsMetricsFile),
                    new ReferenceMetricsChecker(maybeRefSampleConfiguration.get().WgsMetricsFile),
                    new FlagstatChecker(
                            maybeTumorSampleConfiguration.get().FlagstatFile, TUM_PROPORTION_MAPPED, TUM_PROPORTION_DUPLICATE),
                    new FlagstatChecker(
                            maybeRefSampleConfiguration.get().FlagstatFile, REF_PROPORTION_MAPPED, REF_PROPORTION_DUPLICATE),
                    new PurpleChecker(maybeTumorSampleConfiguration.get().SampleName, mPurpleDir));
        }
        else
        {
            throw new IllegalArgumentException(String.format("Illegal combination of arguments: [%s, %s, %s]",
                    maybeRefSampleConfiguration,
                    maybeTumorSampleConfiguration,
                    mPurpleDir));
        }

        List<QCValue> qcValues = Lists.newArrayList();
        for(HealthChecker checker : checkers)
        {
            qcValues.addAll(checker.run());
        }

        for(QCValue qcValue : qcValues)
        {
            HC_LOGGER.info("QC Metric '{}' has value '{}'", qcValue.Type, qcValue.Value);
        }

        if(HealthCheckEvaluation.isPass(qcValues))
        {
            if(mWriteEvaluationFile)
            {
                String evaluationFile = fileOutputBasePath() + ".HealthCheckSucceeded";
                new FileOutputStream(evaluationFile).close();
                HC_LOGGER.info("Evaluation file written to {}", evaluationFile);
            }
            HC_LOGGER.info("Health check evaluation succeeded.");
        }
        else
        {
            if(mWriteEvaluationFile)
            {
                String evaluationFile = fileOutputBasePath() + ".HealthCheckFailed";
                new FileOutputStream(evaluationFile).close();
                HC_LOGGER.info("Evaluation file written to {}", evaluationFile);
            }
            HC_LOGGER.info("Health check evaluation failed!");
        }
    }

    private String fileOutputBasePath()
    {
        String sampleName = mTumorSampleConfig != null ? mTumorSampleConfig.SampleName : mRefSampleConfig.SampleName;
        return mOutputDir + sampleName;
    }

    public static void main(final String... args) throws IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder("HealthChecker");
        registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        HealthChecksApplication healthChecksApplication = new HealthChecksApplication(configBuilder);
        healthChecksApplication.run();
    }

    private static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(REF_SAMPLE, "The name of the reference sample");
        configBuilder.addConfigItem(TUMOR_SAMPLE, "The name of the tumor sample");
        configBuilder.addPath(PURPLE_DIR_CFG, false, PURPLE_DIR_DESC);
        configBuilder.addPath(REF_WGS_METRICS_FILE, false, "The path to the wgs metrics file of reference sample");
        configBuilder.addPath(TUMOR_WGS_METRICS_FILE, false, "The path to the wgs metrics file of tumor sample");
        configBuilder.addPath(REF_FLAGSTAT_FILE, false, "The path to the flagstat file of reference sample");
        configBuilder.addPath(TUMOR_FLAGSTAT_FILE, false, "The path to the flagstat file of tumor sample");
        configBuilder.addFlag(DO_NOT_WRITE_EVALUATION_FILE, "Do not write final success or failure file");
        addOutputDir(configBuilder);
        addLoggingOptions(configBuilder);
    }
}
