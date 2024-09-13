package com.hartwig.hmftools.healthchecker;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REF_METRICS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REF_METRICS_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_METRICS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_METRICS_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import java.io.FileOutputStream;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class HealthChecksApplication
{
    public static final Logger HC_LOGGER = LogManager.getLogger(HealthChecksApplication.class);

    private static final String DO_NOT_WRITE_EVALUATION_FILE = "do_not_write_evaluation_file";

    private final String mTumorId;
    private final String mReferenceId;
    private final String mTumorMetricsDir;
    private final String mReferenceMetricsDir;
    private final String mPurpleDir;
    private final String mOutputDir;

    private final boolean mWriteEvaluationFile;

    public HealthChecksApplication(final ConfigBuilder configBuilder)
    {
        mReferenceId = configBuilder.getValue(REFERENCE);
        mTumorId = configBuilder.getValue(TUMOR);
        mReferenceMetricsDir = configBuilder.getValue(REF_METRICS_DIR_CFG);
        mTumorMetricsDir = configBuilder.getValue(TUMOR_METRICS_DIR_CFG);

        mPurpleDir = configBuilder.getValue(PURPLE_DIR_CFG);
        mOutputDir = parseOutputDir(configBuilder);
        mWriteEvaluationFile = !configBuilder.hasFlag(DO_NOT_WRITE_EVALUATION_FILE);
    }

    public void run()
    {
        List<QCValue> qcValues = Lists.newArrayList();

        if(mTumorId != null && mTumorMetricsDir != null)
            qcValues.addAll(MetricsDataLoader.loadMetricValues(mTumorId, mTumorMetricsDir, true));

        if(mReferenceId != null && mReferenceMetricsDir != null)
            qcValues.addAll(MetricsDataLoader.loadMetricValues(mReferenceId, mReferenceMetricsDir, false));

        String purpleSampleId = mTumorId != null ? mTumorId : mReferenceId;

        qcValues.addAll(PurpleDataLoader.loadPurpleValues(purpleSampleId, mPurpleDir));

        for(QCValue qcValue : qcValues)
        {
            HC_LOGGER.info("QC Metric '{}' has value '{}'", qcValue.Type, qcValue.Value);
        }

        String evaluationFile = fileOutputBasePath();

        if(HealthCheckEvaluation.isPass(qcValues))
        {
            HC_LOGGER.info("Health check evaluation succeeded");
            evaluationFile += ".HealthCheckSucceeded";
        }
        else
        {
            HC_LOGGER.info("Health check evaluation failed");
            evaluationFile += ".HealthCheckFailed";
        }

        if(mWriteEvaluationFile)
        {
            try
            {
                new FileOutputStream(evaluationFile).close();
                HC_LOGGER.info("evaluation file written to {}", evaluationFile);
            }
            catch(IOException e)
            {
                HC_LOGGER.error("failed to write health-checker output: {}", e.toString());
                System.exit(1);
            }
        }
    }

    private String fileOutputBasePath()
    {
        String sampleName = mTumorId != null ? mTumorId : mReferenceId;
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
        configBuilder.addConfigItem(REFERENCE, REFERENCE_DESC);
        configBuilder.addConfigItem(TUMOR, TUMOR_DESC);
        configBuilder.addPath(PURPLE_DIR_CFG, false, PURPLE_DIR_DESC);
        configBuilder.addPath(REF_METRICS_DIR_CFG, false, REF_METRICS_DIR_DESC);
        configBuilder.addPath(TUMOR_METRICS_DIR_CFG, false, TUMOR_METRICS_DIR_DESC);
        configBuilder.addFlag(DO_NOT_WRITE_EVALUATION_FILE, "Do not write final success or failure file");
        addOutputDir(configBuilder);
        addLoggingOptions(configBuilder);
    }
}
