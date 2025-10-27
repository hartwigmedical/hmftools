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
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REF_METRICS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REF_METRICS_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_IDS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_METRICS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_METRICS_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import static common.QseeConstants.QC_LOGGER;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;

import common.SampleType;

public class PrepConfig
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

    public final String OutputDir;
    public final String OutputId;

    public final int Threads;

    public static final String FLD_TUMOR_ID = "TumorId";
    public static final String FLD_REFERENCE_ID = "ReferenceId";

    public PrepConfig(final ConfigBuilder configBuilder)
    {
        TumorIds = new ArrayList<>();
        ReferenceIds = new ArrayList<>();

        if(configBuilder.hasValue(SAMPLE_ID_FILE))
        {
            loadSampleIdsFile(configBuilder, TumorIds, ReferenceIds);
        }
        else
        {
            TumorIds.addAll(loadSampleIdsString(configBuilder.getValue(TUMOR)));

            if(configBuilder.getValue(REFERENCE) != null)
                ReferenceIds.addAll(loadSampleIdsString(configBuilder.getValue(REFERENCE)));
        }

        AmberDir = configBuilder.getValue(AMBER_DIR_CFG);
        CobaltDir = configBuilder.getValue(COBALT_DIR_CFG);
        EsveeDir = configBuilder.getValue(ESVEE_DIR_CFG);
        PurpleDir = configBuilder.getValue(PURPLE_DIR_CFG);
        ReduxDir = configBuilder.getValue(REDUX_DIR_CFG);
        SageDir = configBuilder.getValue(SAGE_DIR_CFG);

        TumorMetricsDir = configBuilder.getValue(TUMOR_METRICS_DIR_CFG);
        RefMetricsDir = configBuilder.getValue(REF_METRICS_DIR_CFG);

        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);

        Threads = TaskExecutor.parseThreads(configBuilder);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(TUMOR, true, TUMOR_IDS_DESC);
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

        FileWriterUtils.addOutputOptions(configBuilder);

        configBuilder.addConfigItem(THREADS, false, THREADS_DESC, "1");

        ConfigUtils.addLoggingOptions(configBuilder);
    }

    private static List<String> loadSampleIdsString(final String sampleIdsString)
    {
        return List.of(sampleIdsString.split(","));
    }

    private static void loadSampleIdsFile(final ConfigBuilder configBuilder, final List<String> tumorIds, List<String> referenceIds)
    {
        try
        {
            List<String> lines = Files.readAllLines(Paths.get(configBuilder.getValue(SAMPLE_ID_FILE)));

            String header = lines.get(0);
            Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
            Integer tumorIdIndex = fieldsIndexMap.get(FLD_TUMOR_ID);
            Integer referenceIdIndex = fieldsIndexMap.get(FLD_REFERENCE_ID);

            lines.remove(0);

            for(String line : lines)
            {
                String[] values = line.split(TSV_DELIM, -1);
                tumorIds.add(values[tumorIdIndex]);

                if(referenceIdIndex != null)
                    referenceIds.add(values[referenceIdIndex]);
            }

            QC_LOGGER.info("Loaded {} tumor samples and {} reference samples", tumorIds.size(), referenceIds.size());
        }
        catch(IOException e)
        {
            QC_LOGGER.error("Failed to load sample ids file: {}", e.toString());
            e.printStackTrace();
            System.exit(1);
        }
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
