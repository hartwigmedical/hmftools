package prep;

import static com.hartwig.hmftools.common.perf.TaskExecutor.THREADS;
import static com.hartwig.hmftools.common.perf.TaskExecutor.THREADS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.AMBER_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.AMBER_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.COBALT_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.COBALT_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.ESVEE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.ESVEE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REDUX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REDUX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REF_METRICS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REF_METRICS_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_METRICS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_METRICS_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import static common.QSeeConstants.QC_LOGGER;

import java.util.List;

import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;

public class PrepConfig
{
    public final List<String> SampleIds;
    public final String SampleDataDir;

    public final String AmberDir;
    public final String CobaltDir;
    public final String EsveeDir;
    public final String PurpleDir;
    public final String ReduxDir;
    public final String SageDir;
    public final String TumorMetricsDir;
    public final String RefMetricsDir;

    public final String OutputDir;
    public final String OutputId;

    public final int Threads;

    public PrepConfig(final ConfigBuilder configBuilder)
    {
        SampleIds = parseSampleConfig(configBuilder);
        SampleDataDir = configBuilder.getValue(SAMPLE_DATA_DIR_CFG, "");

        AmberDir = configBuilder.getValue(AMBER_DIR_CFG, SampleDataDir);
        CobaltDir = configBuilder.getValue(COBALT_DIR_CFG, SampleDataDir);
        EsveeDir = configBuilder.getValue(ESVEE_DIR_CFG, SampleDataDir);
        PurpleDir = configBuilder.getValue(PURPLE_DIR_CFG, SampleDataDir);
        ReduxDir = configBuilder.getValue(REDUX_DIR_CFG, SampleDataDir);
        SageDir = configBuilder.getValue(SAGE_DIR_CFG, SampleDataDir);

        TumorMetricsDir = configBuilder.getValue(TUMOR_METRICS_DIR_CFG, SampleDataDir);
        RefMetricsDir = configBuilder.getValue(REF_METRICS_DIR_CFG, SampleDataDir);

        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);

        Threads = TaskExecutor.parseThreads(configBuilder);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, false, SAMPLE_DESC);
        configBuilder.addPath(SAMPLE_ID_FILE, false, SAMPLE_ID_FILE_DESC);

        configBuilder.addPath(SAMPLE_DATA_DIR_CFG, false, SAMPLE_DATA_DIR_DESC);

        configBuilder.addPath(AMBER_DIR_CFG, false, AMBER_DIR_DESC);
        configBuilder.addPath(COBALT_DIR_CFG, false, COBALT_DIR_DESC);
        configBuilder.addPath(ESVEE_DIR_CFG, false, ESVEE_DIR_DESC);
        configBuilder.addPath(PURPLE_DIR_CFG, false, PURPLE_DIR_DESC);
        configBuilder.addPath(REDUX_DIR_CFG, false, REDUX_DIR_DESC);
        configBuilder.addPath(SAGE_DIR_CFG, false, SAGE_DIR_DESC);
        configBuilder.addPath(TUMOR_METRICS_DIR_CFG, false, TUMOR_METRICS_DIR_DESC);
        configBuilder.addPath(REF_METRICS_DIR_CFG, false, REF_METRICS_DIR_DESC);

        FileWriterUtils.addOutputOptions(configBuilder);

        configBuilder.addConfigItem(THREADS, false, THREADS_DESC, "1");

        ConfigUtils.addLoggingOptions(configBuilder);
    }

    private static List<String> parseSampleConfig(ConfigBuilder configBuilder)
    {
        if(configBuilder.hasValue(SAMPLE))
            return List.of(configBuilder.getValue(SAMPLE));

        if(configBuilder.hasValue(SAMPLE_ID_FILE))
            return ConfigUtils.loadSampleIdsFile(configBuilder);

        QC_LOGGER.error("Either -{} or -{} must be provided", SAMPLE, SAMPLE_ID_FILE);
        System.exit(1);
        return null;
    }

    public boolean isMultiSample() { return SampleIds.size() > 1; }
    public boolean isSingleSample() { return SampleIds.size() == 1; }

    public String getAmberDir(final String sampleId) { return convertWildcardSamplePath(AmberDir, sampleId); }
    public String getCobaltDir(final String sampleId) { return convertWildcardSamplePath(CobaltDir, sampleId); }
    public String getEsveeDir(final String sampleId) { return convertWildcardSamplePath(EsveeDir, sampleId); }
    public String getPurpleDir(final String sampleId) { return convertWildcardSamplePath(PurpleDir, sampleId); }
    public String getReduxDir(final String sampleId) { return convertWildcardSamplePath(ReduxDir, sampleId); }
    public String getSageDir(final String sampleId) { return convertWildcardSamplePath(SageDir, sampleId); }
    public String getTumorMetricsDir(final String sampleId) { return convertWildcardSamplePath(TumorMetricsDir, sampleId); }
    public String getRefMetricsDir(final String sampleId) { return convertWildcardSamplePath(RefMetricsDir, sampleId); }
}
