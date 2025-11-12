package com.hartwig.hmftools.qsee.prep;

import static com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL;
import static com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_DESC;
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
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REF_METRICS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REF_METRICS_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_METRICS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_METRICS_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import java.util.List;

import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

import com.hartwig.hmftools.qsee.common.SampleIdsLoader;
import com.hartwig.hmftools.qsee.common.SampleType;

public class CommonPrepConfig
{
    public final List<String> TumorIds;
    public final List<String> ReferenceIds;

    public final String AmberDir;
    public final String CobaltDir;
    public final String EsveeDir;
    public final String PurpleDir;
    public final String ReduxDir;
    public final String SageDir;
    public final String TumorMetricsDir;
    public final String RefMetricsDir;

    public final List<DriverGene> DriverGenes;

    public final boolean AllowMissingInput;

    public final String OutputDir;

    public final int Threads;

    public static final String ALLOW_MISSING_INPUT = "allow_missing_input";
    public static final String ALLOW_MISSING_INPUT_DESC = "Continue sample data extraction even if some input files are missing";

    public CommonPrepConfig(final ConfigBuilder configBuilder)
    {
        SampleIdsLoader sampleIdsLoader = new SampleIdsLoader().fromConfig(configBuilder);
        TumorIds = sampleIdsLoader.tumorIds();
        ReferenceIds = sampleIdsLoader.referenceIds();

        AmberDir = configBuilder.getValue(AMBER_DIR_CFG);
        CobaltDir = configBuilder.getValue(COBALT_DIR_CFG);
        EsveeDir = configBuilder.getValue(ESVEE_DIR_CFG);
        PurpleDir = configBuilder.getValue(PURPLE_DIR_CFG);
        ReduxDir = configBuilder.getValue(REDUX_DIR_CFG);
        SageDir = configBuilder.getValue(SAGE_DIR_CFG);
        TumorMetricsDir = configBuilder.getValue(TUMOR_METRICS_DIR_CFG);
        RefMetricsDir = configBuilder.getValue(REF_METRICS_DIR_CFG);

        DriverGenes = DriverGenePanelConfig.loadDriverGenes(configBuilder);

        AllowMissingInput = configBuilder.hasFlag(ALLOW_MISSING_INPUT);

        OutputDir = parseOutputDir(configBuilder);

        Threads = TaskExecutor.parseThreads(configBuilder);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        // TODO: Set required to true for some args

        configBuilder.addConfigItem(TUMOR, false, TUMOR_IDS_DESC);
        configBuilder.addConfigItem(REFERENCE, false, REFERENCE_IDS_DESC);
        configBuilder.addPath(SAMPLE_ID_FILE, false, SAMPLE_ID_FILE_DESC);

        configBuilder.addPath(AMBER_DIR_CFG, false, AMBER_DIR_DESC);
        configBuilder.addPath(COBALT_DIR_CFG, false, COBALT_DIR_DESC);
        configBuilder.addPath(ESVEE_DIR_CFG, false, ESVEE_DIR_DESC);
        configBuilder.addPath(PURPLE_DIR_CFG, false, PURPLE_DIR_DESC);
        configBuilder.addPath(REDUX_DIR_CFG, false, REDUX_DIR_DESC);
        configBuilder.addPath(SAGE_DIR_CFG, false, SAGE_DIR_DESC);
        configBuilder.addPath(TUMOR_METRICS_DIR_CFG, false, TUMOR_METRICS_DIR_DESC);
        configBuilder.addPath(REF_METRICS_DIR_CFG, false, REF_METRICS_DIR_DESC);

        configBuilder.addPath(DRIVER_GENE_PANEL, false, DRIVER_GENE_PANEL_DESC);

        configBuilder.addFlag(ALLOW_MISSING_INPUT, ALLOW_MISSING_INPUT_DESC);

        configBuilder.addPath(OUTPUT_DIR, true, OUTPUT_DIR_DESC);

        configBuilder.addConfigItem(THREADS, false, THREADS_DESC, "1");
        ConfigUtils.addLoggingOptions(configBuilder);
    }

    public List<String> getSampleIds(SampleType sampleType)
    {
        return sampleType == SampleType.TUMOR ? TumorIds : ReferenceIds;
    }

    public String getAmberDir(final String sampleId) { return convertWildcardSamplePath(AmberDir, sampleId); }
    public String getCobaltDir(final String sampleId) { return convertWildcardSamplePath(CobaltDir, sampleId); }
    public String getEsveeDir(final String sampleId) { return convertWildcardSamplePath(EsveeDir, sampleId); }
    public String getPurpleDir(final String sampleId) { return convertWildcardSamplePath(PurpleDir, sampleId); }
    public String getReduxDir(final String sampleId) { return convertWildcardSamplePath(ReduxDir, sampleId); }
    public String getSageDir(final String sampleId) { return convertWildcardSamplePath(SageDir, sampleId); }

    public String getBamMetricsDir(final String sampleId, final SampleType sampleType)
    {
        String baseDir = sampleType == SampleType.TUMOR ? TumorMetricsDir : RefMetricsDir;
        return convertWildcardSamplePath(baseDir, sampleId);
    }
}
