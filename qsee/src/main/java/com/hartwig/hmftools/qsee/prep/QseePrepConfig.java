package com.hartwig.hmftools.qsee.prep;

import static com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL;
import static com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_DESC;
import static com.hartwig.hmftools.common.perf.TaskExecutor.THREADS;
import static com.hartwig.hmftools.common.perf.TaskExecutor.THREADS_DESC;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SEQUENCING_TYPE_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.COBALT_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.COBALT_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.ESVEE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.ESVEE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REDUX_REF_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REDUX_REF_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REDUX_TUMOR_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REDUX_TUMOR_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REF_METRICS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REF_METRICS_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_METRICS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_METRICS_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.qsee.cohort.CohortPercentilesFile.COHORT_PERCENTILES_FILE_CFG;
import static com.hartwig.hmftools.qsee.cohort.CohortPercentilesFile.COHORT_PERCENTILES_FILE_CFG_DESC;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;
import static com.hartwig.hmftools.qsee.status.ThresholdOverridesFile.THRESHOLD_OVERRIDES_FILE_CFG;
import static com.hartwig.hmftools.qsee.status.ThresholdOverridesFile.THRESHOLD_OVERRIDES_FILE_CFG_DESC;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

import com.hartwig.hmftools.qsee.common.SampleIdsLoader;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.prep.category.PrepCategory;
import com.hartwig.hmftools.qsee.status.ThresholdOverridesFile;
import com.hartwig.hmftools.qsee.status.ThresholdRegistry;

public class QseePrepConfig
{
    public final List<String> TumorIds;
    public final List<String> ReferenceIds;

    public final String SampleDataDir;

    public final String BamMetricsDirTumor;
    public final String BamMetricsDirRef;
    public final String ReduxDirTumor;
    public final String ReduxDirRef;

    public final String CobaltDir;
    public final String EsveeDir;
    public final String PurpleDir;
    public final String SageDir;

    public final List<DriverGene> DriverGenes;
    public final String CohortPercentilesFile;
    public final String ThresholdsFile;
    public final ThresholdRegistry QcThresholds;

    public final SequencingType SequencingTech;
    public final List<PrepCategory> Categories;
    public final boolean AllowMissingInput;

    public final String OutputDir;
    public final String OutputId;

    public final int Threads;

    public static final String SAMPLE_ID_FILE_DESC = "Sample ID TSV file containing columns TumorId and optionally ReferenceId";

    public static final String SKIP_CATEGORIES = "skip_categories";
    public static final String SKIP_CATEGORIES_DESC = "Comma-separated list of categories to skip";
    public static final String SKIP_CATEGORIES_DELIM = ",";

    public static final String ALLOW_MISSING_INPUT = "allow_missing_input";
    public static final String ALLOW_MISSING_INPUT_DESC = "Continue sample data extraction even if some input files are missing";

    public QseePrepConfig(final ConfigBuilder configBuilder)
    {
        SampleIdsLoader sampleIdsLoader = new SampleIdsLoader().fromConfig(configBuilder);
        TumorIds = sampleIdsLoader.tumorIds();
        ReferenceIds = sampleIdsLoader.referenceIds();

        SampleDataDir = configBuilder.getValue(SAMPLE_DATA_DIR_CFG, "");

        BamMetricsDirTumor = configBuilder.getValue(TUMOR_METRICS_DIR_CFG, SampleDataDir);
        BamMetricsDirRef = configBuilder.getValue(REF_METRICS_DIR_CFG, SampleDataDir);
        ReduxDirTumor = configBuilder.getValue(REDUX_TUMOR_DIR_CFG, SampleDataDir);
        ReduxDirRef = configBuilder.getValue(REDUX_REF_DIR_CFG, SampleDataDir);

        CobaltDir = configBuilder.getValue(COBALT_DIR_CFG, SampleDataDir);
        EsveeDir = configBuilder.getValue(ESVEE_DIR_CFG, SampleDataDir);
        PurpleDir = configBuilder.getValue(PURPLE_DIR_CFG, SampleDataDir);
        SageDir = configBuilder.getValue(SAGE_DIR_CFG, SampleDataDir);

        DriverGenes = DriverGenePanelConfig.loadDriverGenes(configBuilder);
        CohortPercentilesFile = configBuilder.getValue(COHORT_PERCENTILES_FILE_CFG);
        ThresholdsFile = configBuilder.getValue(THRESHOLD_OVERRIDES_FILE_CFG);

        QcThresholds = ThresholdsFile == null
                ? ThresholdRegistry.createDefault()
                : ThresholdOverridesFile.read(ThresholdsFile);

        SequencingTech = SequencingType.valueOf(configBuilder.getValue(SEQUENCING_TYPE_CFG));
        Categories = parseCategories(configBuilder.getValue(SKIP_CATEGORIES));
        AllowMissingInput = configBuilder.hasFlag(ALLOW_MISSING_INPUT);

        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);

        Threads = TaskExecutor.parseThreads(configBuilder);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(TUMOR, false, TUMOR_IDS_DESC);
        configBuilder.addConfigItem(REFERENCE, false, REFERENCE_IDS_DESC);
        configBuilder.addPath(SAMPLE_ID_FILE, false, SAMPLE_ID_FILE_DESC);

        configBuilder.addPath(SAMPLE_DATA_DIR_CFG, false, SAMPLE_DATA_DIR_DESC);

        configBuilder.addPath(TUMOR_METRICS_DIR_CFG, false, TUMOR_METRICS_DIR_DESC);
        configBuilder.addPath(REF_METRICS_DIR_CFG, false, REF_METRICS_DIR_DESC);
        configBuilder.addPath(REDUX_TUMOR_DIR_CFG, false, REDUX_TUMOR_DIR_DESC);
        configBuilder.addPath(REDUX_REF_DIR_CFG, false, REDUX_REF_DIR_DESC);

        configBuilder.addPath(COBALT_DIR_CFG, false, COBALT_DIR_DESC);
        configBuilder.addPath(ESVEE_DIR_CFG, false, ESVEE_DIR_DESC);
        configBuilder.addPath(PURPLE_DIR_CFG, false, PURPLE_DIR_DESC);
        configBuilder.addPath(SAGE_DIR_CFG, false, SAGE_DIR_DESC);

        SequencingType.registerConfig(configBuilder);
        configBuilder.addPath(DRIVER_GENE_PANEL, false, DRIVER_GENE_PANEL_DESC);
        configBuilder.addPath(COHORT_PERCENTILES_FILE_CFG, false, COHORT_PERCENTILES_FILE_CFG_DESC);
        configBuilder.addPath(THRESHOLD_OVERRIDES_FILE_CFG, false, THRESHOLD_OVERRIDES_FILE_CFG_DESC);

        configBuilder.addConfigItem(SEQUENCING_TYPE_CFG, false, SequencingType.ILLUMINA.toString());
        configBuilder.addConfigItem(SKIP_CATEGORIES, false, SKIP_CATEGORIES_DESC, null);
        configBuilder.addFlag(ALLOW_MISSING_INPUT, ALLOW_MISSING_INPUT_DESC);

        configBuilder.addPath(OUTPUT_DIR, true, OUTPUT_DIR_DESC);
        configBuilder.addConfigItem(OUTPUT_ID, false, OUTPUT_DIR_DESC);

        configBuilder.addConfigItem(THREADS, false, THREADS_DESC, "1");
        ConfigUtils.addLoggingOptions(configBuilder);
    }

    public List<String> getSampleIds(SampleType sampleType) { return sampleType == SampleType.TUMOR ? TumorIds : ReferenceIds; }

    public boolean isSinglePatient() { return TumorIds.size() <= 1 && ReferenceIds.size() <= 1; }
    public String getCobaltDir(String sampleId) { return convertWildcardSamplePath(CobaltDir, sampleId); }
    public String getEsveeDir(String sampleId) { return convertWildcardSamplePath(EsveeDir, sampleId); }
    public String getPurpleDir(String sampleId) { return convertWildcardSamplePath(PurpleDir, sampleId); }
    public String getSageDir(String sampleId) { return convertWildcardSamplePath(SageDir, sampleId); }

    public String getBamMetricsDir(String sampleId, SampleType sampleType)
    {
        String baseDir = sampleType == SampleType.TUMOR ? BamMetricsDirTumor : BamMetricsDirRef;
        return convertWildcardSamplePath(baseDir, sampleId);
    }

    public String getReduxDir(String sampleId, SampleType sampleType)
    {
        String baseDir = sampleType == SampleType.TUMOR ? ReduxDirTumor : ReduxDirRef;
        return convertWildcardSamplePath(baseDir, sampleId);
    }

    public static List<PrepCategory> parseCategories(String skipCategoriesString)
    {
        if(skipCategoriesString == null)
            return List.of(PrepCategory.values());

        List<PrepCategory> requiredCategories = List.of(PrepCategory.SUMMARY_TABLE_PURPLE, PrepCategory.SUMMARY_TABLE_BAM_METRICS);
        List<PrepCategory> categoriesToSkip = new ArrayList<>();
        for(String string : skipCategoriesString.split(SKIP_CATEGORIES_DELIM))
        {
            PrepCategory categoryToSkip = PrepCategory.valueOf(string.toUpperCase());

            if(requiredCategories.contains(categoryToSkip))
            {
                QC_LOGGER.error("Skipping category({}) is not supported", categoryToSkip);
                System.exit(1);
            }

            categoriesToSkip.add(categoryToSkip);
        }

        return Stream.of(PrepCategory.values()).filter(c -> !categoriesToSkip.contains(c)).toList();
    }
}
