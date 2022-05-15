package com.hartwig.hmftools.isofox.loader;

import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadGeneIdsFile;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.isofox.IsofoxConfig.GENE_ID_FILE;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.ITEM_DELIM;
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
import com.hartwig.hmftools.common.utils.FileReaderUtils;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class DataLoaderConfig
{
    public static final String SAMPLE = "sample";
    public static final String SAMPLE_ID_FILE = "sample_id_file";

    public static final String SAMPLE_DATA_DIRECTORY = "sample_data_dir";
    public static final String GENE_DATA_DIRECTORY = "gene_data_dir";
    public static final String ALT_SJ_DATA_DIRECTORY = "alt_sj_data_dir";
    public static final String FUSION_DATA_DIRECTORY = "fusion_data_dir";
    public static final String STATISTICS_DATA_DIRECTORY = "stats_data_dir";

    public static final String CANCER_TYPES_FILE = "cancer_types_file";
    public static final String GENE_DIST_FILE = "gene_distribution_file";
    public static final String ALT_SJ_COHORT_FILE = "alt_sj_cohort_file";

    public static final String FLD_SAMPLE_ID = "SampleId";
    public static final String FLD_CANCER_TYPE = "CancerType";

    public static final String LOAD_TYPES = "load_types";
    public static final String THREADS = "threads";

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

    public final CommandLine CmdLineArgs;

    public DataLoaderConfig(final CommandLine cmd)
    {
        CmdLineArgs = cmd;
        SampleIds = Lists.newArrayList();
        PrimaryCancerTypes = Lists.newArrayList();
        LoadTypes = Lists.newArrayList();
        SampleCancerTypes = Maps.newHashMap();

        if(cmd.hasOption(SAMPLE))
        {
            addSampleInfo(cmd.getOptionValue(SAMPLE));
        }
        else if(cmd.hasOption(SAMPLE_ID_FILE))
        {
            loadSampleDataFile(cmd.getOptionValue(SAMPLE_ID_FILE));
        }

        if(cmd.hasOption(LOAD_TYPES))
        {
            LoadTypes.addAll(Arrays.stream(cmd.getOptionValue(LOAD_TYPES)
                    .split(ITEM_DELIM, -1))
                    .map(x -> DataLoadType.valueOf(x))
                    .collect(Collectors.toList()));

            ISF_LOGGER.info("loading types: {}", LoadTypes);
        }

        if(cmd.hasOption(CANCER_TYPES_FILE))
        {
            try
            {
                final List<String> lines = Files.readAllLines(Paths.get(cmd.getOptionValue(CANCER_TYPES_FILE)));
                lines.stream().filter(x -> !x.equals("CancerType")).forEach(x -> PrimaryCancerTypes.add(x));
                ISF_LOGGER.info("loaded {} known cancer types", PrimaryCancerTypes.size());
            }
            catch(IOException e)
            {
                ISF_LOGGER.warn("invalid cancer types file: {}", e.toString());
            }
        }

        SampleDataDir = checkAddDirSeparator(cmd.getOptionValue(SAMPLE_DATA_DIRECTORY));
        GeneDistributionFile = cmd.getOptionValue(GENE_DIST_FILE);
        AltSjCohortFile = cmd.getOptionValue(ALT_SJ_COHORT_FILE);

        GeneDataDir = cmd.hasOption(GENE_DATA_DIRECTORY) ?
                checkAddDirSeparator(cmd.getOptionValue(GENE_DATA_DIRECTORY)) : SampleDataDir;

        AltSjDataDir = cmd.hasOption(ALT_SJ_DATA_DIRECTORY) ?
                checkAddDirSeparator(cmd.getOptionValue(ALT_SJ_DATA_DIRECTORY)) : SampleDataDir;

        FusionDataDir = cmd.hasOption(FUSION_DATA_DIRECTORY) ?
                checkAddDirSeparator(cmd.getOptionValue(FUSION_DATA_DIRECTORY)) : SampleDataDir;

        StatisticsDataDir = cmd.hasOption(STATISTICS_DATA_DIRECTORY) ?
                checkAddDirSeparator(cmd.getOptionValue(STATISTICS_DATA_DIRECTORY)) : SampleDataDir;

        RestrictedGeneIds = Lists.newArrayList();

        Threads = Integer.parseInt(cmd.getOptionValue(THREADS, "0"));

        if(cmd.hasOption(GENE_ID_FILE))
        {
            final String inputFile = cmd.getOptionValue(GENE_ID_FILE);
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

    public static Options createCmdLineOptions()
    {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "Sample tumor ID");
        options.addOption(SAMPLE_ID_FILE, true, "File with list of samples and cancer types to load data for");
        options.addOption(LOAD_TYPES, true, "Load specific types only (default=ALL)");
        options.addOption(SAMPLE_DATA_DIRECTORY, true, "Sample data directory");
        options.addOption(GENE_DATA_DIRECTORY, true, "Gene data directory, will use sample data dir if not present");
        options.addOption(ALT_SJ_DATA_DIRECTORY, true, "Alt-SJ data directory, will use sample data dir if not present");
        options.addOption(FUSION_DATA_DIRECTORY, true, "Fusion data directory, will use sample data dir if not present");
        options.addOption(STATISTICS_DATA_DIRECTORY, true, "Summary statistics data directory, will use sample data dir if not present");
        options.addOption(CANCER_TYPES_FILE, true, "Primary cancer types (otherwise will use 'Other' for sample");
        options.addOption(GENE_DIST_FILE, true, "Gene distribution for medians and percentile data");
        options.addOption(ALT_SJ_COHORT_FILE, true, "Alternate splice junction cohort file");
        options.addOption(GENE_ID_FILE, true, "Optional CSV file of genes to analyse");
        options.addOption(THREADS, true, "Optional multi-threading when loading multiple samples");
        options.addOption(LOG_DEBUG, false, "Log verbose");

        addDatabaseCmdLineArgs(options);

        return options;
    }
}
