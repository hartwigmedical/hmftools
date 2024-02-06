package com.hartwig.hmftools.isofox.loader;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.GENE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.GENE_ID_FILE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.IGNORE_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadGeneIdsFile;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CANCER_TYPE;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.file.FileReaderUtils;

public class DataLoaderConfig
{
    public static final String GENE_DATA_DIRECTORY = "gene_data_dir";
    public static final String ALT_SJ_DATA_DIRECTORY = "alt_sj_data_dir";
    public static final String FUSION_DATA_DIRECTORY = "fusion_data_dir";
    public static final String STATISTICS_DATA_DIRECTORY = "stats_data_dir";

    public static final String CANCER_TYPES_FILE = "cancer_types_file";
    public static final String GENE_DIST_FILE = "gene_distribution_file";
    public static final String GENE_DIST_FILE_DESC = "Gene distribution for medians and percentile data";
    public static final String ALT_SJ_COHORT_FILE = "alt_sj_cohort_file";

    public static final String LOAD_TYPES = "load_types";

    public final String SampleDataDir;
    public final String GeneDataDir;
    public final String AltSjDataDir;
    public final String FusionDataDir;
    public final String StatisticsDataDir;
    public final String GeneDistributionFile;
    public final String AltSjCohortFile;
    public final List<String> RestrictedGeneIds;
    public final List<String> SampleIds;
    public final List<String> PrimaryCancerTypes;
    public final List<DataLoadType> LoadTypes;
    public final Map<String,String> SampleCancerTypes;
    public final Integer Threads;

    public final ConfigBuilder ConfigItems;

    public DataLoaderConfig(final ConfigBuilder configBuilder)
    {
        ConfigItems = configBuilder;
        SampleIds = Lists.newArrayList();
        PrimaryCancerTypes = Lists.newArrayList();
        LoadTypes = Lists.newArrayList();
        SampleCancerTypes = Maps.newHashMap();

        if(configBuilder.hasValue(SAMPLE))
        {
            addSampleInfo(configBuilder.getValue(SAMPLE));
        }
        else if(configBuilder.hasValue(SAMPLE_ID_FILE))
        {
            loadSampleDataFile(configBuilder.getValue(SAMPLE_ID_FILE));
        }

        if(configBuilder.hasValue(LOAD_TYPES))
        {
            LoadTypes.addAll(Arrays.stream(configBuilder.getValue(LOAD_TYPES)
                    .split(ITEM_DELIM, -1))
                    .map(x -> DataLoadType.valueOf(x))
                    .collect(Collectors.toList()));

            ISF_LOGGER.info("loading types: {}", LoadTypes);
        }

        if(configBuilder.hasValue(CANCER_TYPES_FILE))
        {
            try
            {
                final List<String> lines = Files.readAllLines(Paths.get(configBuilder.getValue(CANCER_TYPES_FILE)));
                lines.stream().filter(x -> !x.equals("CancerType")).forEach(x -> PrimaryCancerTypes.add(x));
                ISF_LOGGER.info("loaded {} known cancer types", PrimaryCancerTypes.size());
            }
            catch(IOException e)
            {
                ISF_LOGGER.warn("invalid cancer types file: {}", e.toString());
            }
        }

        SampleDataDir = checkAddDirSeparator(configBuilder.getValue(SAMPLE_DATA_DIR_CFG));
        GeneDistributionFile = configBuilder.getValue(GENE_DIST_FILE);
        AltSjCohortFile = configBuilder.getValue(ALT_SJ_COHORT_FILE);

        GeneDataDir = configBuilder.hasValue(GENE_DATA_DIRECTORY) ?
                checkAddDirSeparator(configBuilder.getValue(GENE_DATA_DIRECTORY)) : SampleDataDir;

        AltSjDataDir = configBuilder.hasValue(ALT_SJ_DATA_DIRECTORY) ?
                checkAddDirSeparator(configBuilder.getValue(ALT_SJ_DATA_DIRECTORY)) : SampleDataDir;

        FusionDataDir = configBuilder.hasValue(FUSION_DATA_DIRECTORY) ?
                checkAddDirSeparator(configBuilder.getValue(FUSION_DATA_DIRECTORY)) : SampleDataDir;

        StatisticsDataDir = configBuilder.hasValue(STATISTICS_DATA_DIRECTORY) ?
                checkAddDirSeparator(configBuilder.getValue(STATISTICS_DATA_DIRECTORY)) : SampleDataDir;

        RestrictedGeneIds = Lists.newArrayList();

        Threads = parseThreads(configBuilder);

        if(configBuilder.hasValue(GENE_ID_FILE))
        {
            final String inputFile = configBuilder.getValue(GENE_ID_FILE);
            RestrictedGeneIds.addAll(loadGeneIdsFile(inputFile));
            ISF_LOGGER.info("file({}) loaded {} restricted genes", inputFile, RestrictedGeneIds.size());
        }
    }

    public boolean loadDataType(final DataLoadType type)
    {
        return LoadTypes.isEmpty() || LoadTypes.contains(type);
    }

    private void loadSampleDataFile(final String filename)
    {
        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));

            Map<String,Integer> fieldsIndexMap = FileReaderUtils.createFieldsIndexMap(lines.get(0), DELIMITER);
            lines.remove(0);

            if(!fieldsIndexMap.containsKey(FLD_SAMPLE_ID))
            {
                ISF_LOGGER.error("sample ID file missing 'SampleId'");
                return;
            }

            int sampleIdIndex = fieldsIndexMap.get(FLD_SAMPLE_ID);
            Integer cancerTypeIndex = fieldsIndexMap.get(FLD_CANCER_TYPE);

            for(String line : lines)
            {
                if(line.isEmpty() || line.startsWith(IGNORE_SAMPLE_ID))
                    continue;

                String[] values = line.split(DELIMITER, -1);
                String sampleId = values[sampleIdIndex];
                SampleIds.add(sampleId);

                if(cancerTypeIndex != null)
                {
                    String cancerType = values[cancerTypeIndex];
                    SampleCancerTypes.put(sampleId, cancerType);
                }
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.warn("invalid sampleId file: {}", e.toString());
        }
    }

    private void addSampleInfo(final String sampleInfo)
    {
        if(!sampleInfo.contains(DELIMITER))
        {
            SampleIds.add(sampleInfo);
            return;
        }

        String[] items = sampleInfo.split(DELIMITER, -1);
        String sampleId = items[0];
        String cancerType = items[1];
        SampleIds.add(sampleId);
        SampleCancerTypes.put(sampleId, cancerType);
    }

    public boolean processGeneId(final String geneId)
    {
        return RestrictedGeneIds.isEmpty() || RestrictedGeneIds.contains(geneId);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, SAMPLE_DESC);
        configBuilder.addConfigItem(SAMPLE_ID_FILE, "File with list of samples and cancer types to load data for");
        configBuilder.addConfigItem(LOAD_TYPES, "Load specific types only (default=ALL)");
        configBuilder.addPath(SAMPLE_DATA_DIR_CFG, false, SAMPLE_DATA_DIR_DESC);
        configBuilder.addPath(GENE_DATA_DIRECTORY, false, "Gene data directory, will use sample data dir if not present");
        configBuilder.addPath(ALT_SJ_DATA_DIRECTORY, false, "Alt-SJ data directory, will use sample data dir if not present");
        configBuilder.addPath(FUSION_DATA_DIRECTORY, false, "Fusion data directory, will use sample data dir if not present");
        configBuilder.addPath(STATISTICS_DATA_DIRECTORY, false, "Summary statistics data directory, will use sample data dir if not present");
        configBuilder.addPath(CANCER_TYPES_FILE, false, "Primary cancer types (otherwise will use 'Other' for sample");
        configBuilder.addPath(GENE_DIST_FILE, false, GENE_DIST_FILE_DESC);
        configBuilder.addPath(ALT_SJ_COHORT_FILE, false, "Alternate splice junction cohort file");
        configBuilder.addPath(GENE_ID_FILE, false, GENE_ID_FILE_DESC);

        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
        addDatabaseCmdLineArgs(configBuilder, true);
    }
}
