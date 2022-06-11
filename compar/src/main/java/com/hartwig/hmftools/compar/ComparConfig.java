package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION_DESC;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.compar.Category.ALL_CATEGORIES;
import static com.hartwig.hmftools.compar.Category.DRIVER;
import static com.hartwig.hmftools.compar.Category.FUSION;
import static com.hartwig.hmftools.compar.CommonUtils.DATA_DELIM;
import static com.hartwig.hmftools.compar.CommonUtils.ITEM_DELIM;
import static com.hartwig.hmftools.compar.CommonUtils.SUB_ITEM_DELIM;
import static com.hartwig.hmftools.compar.FileSources.fromConfig;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.DB_DEFAULT_ARGS;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;

import java.io.IOException;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.ensemblcache.GeneNameMapping;
import com.hartwig.hmftools.common.utils.ConfigUtils;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ComparConfig
{
    public final List<String> SampleIds;

    public final Map<Category,MatchLevel> Categories;

    public final List<String> SourceNames; // list of sources to compare, eg prod vs pilot, or pipeline_1 vs pipeline_2

    public final Map<String,DatabaseAccess> DbConnections; // database access details keyed by source
    public final Map<String,FileSources> FileSources; // directories per type and keyed by source
    public final Map<String,String> SourceSampleIds; // as required, mapping from sampleId to source sampleId

    public final Set<String> DriverGenes;
    public final GeneNameMapping GeneMapping;

    public final DiffThresholds Thresholds;

    public final String OutputDir;
    public final String OutputId;

    public final boolean WriteDetailed;
    public final int Threads;

    private boolean mIsValid;

    // config strings
    public static final String CATEGORIES = "categories";
    public static final String MATCH_LEVEL = "match_level";

    public static final String DB_SOURCES = "db_sources";
    public static final String FILE_SOURCES = "file_sources";
    public static final String THRESHOLDS = "thresholds";

    public static final String SAMPLE = "sample";
    public static final String SOURCE_SAMPLE_MAPPINGS = "source_sample_mappings";
    public static final String SAMPLE_ID_FILE = "sample_id_file";
    public static final String THREADS = "threads";
    public static final String WRITE_DETAILED_FILES = "write_detailed";

    public static final Logger CMP_LOGGER = LogManager.getLogger(ComparConfig.class);

    public ComparConfig(final CommandLine cmd)
    {
        mIsValid = true;

        SampleIds = Lists.newArrayList();
        loadSampleIds(cmd);

        Categories = Maps.newHashMap();

        MatchLevel matchLevel = MatchLevel.valueOf(cmd.getOptionValue(MATCH_LEVEL, REPORTABLE.toString()));

        if(!cmd.hasOption(CATEGORIES) || cmd.getOptionValue(CATEGORIES).equals(ALL_CATEGORIES))
        {
            Arrays.stream(Category.values()).forEach(x -> Categories.put(x, matchLevel));
        }
        else
        {
            final String[] catDataList = cmd.getOptionValue(CATEGORIES).split(DATA_DELIM);

            for(String catData : catDataList)
            {
                Category category;
                MatchLevel specificMatchLevel;

                if(catData.contains("="))
                {
                    String[] catItems = catData.split("=");
                    category = Category.valueOf(catItems[0]);
                    specificMatchLevel = MatchLevel.valueOf(catItems[1]);
                }
                else
                {
                    category = Category.valueOf(catData);
                    specificMatchLevel = matchLevel;
                }

                Categories.put(category, specificMatchLevel);
            }
        }

        CMP_LOGGER.info("comparing categories: {}", Categories.isEmpty() ? ALL_CATEGORIES : Categories.toString());

        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID);
        WriteDetailed = cmd.hasOption(WRITE_DETAILED_FILES);
        Threads = Integer.parseInt(cmd.getOptionValue(THREADS, "0"));

        DbConnections = Maps.newHashMap();
        FileSources = Maps.newHashMap();
        SourceNames = Lists.newArrayList();
        loadDatabaseSources(cmd);
        loadFileSources(cmd);

        Thresholds = new DiffThresholds();
        Thresholds.loadConfig(cmd.getOptionValue(THRESHOLDS, ""));

        SourceSampleIds = Maps.newHashMap();
        if(cmd.hasOption(SOURCE_SAMPLE_MAPPINGS))
        {
            String[] sampleMappings = cmd.getOptionValue(SOURCE_SAMPLE_MAPPINGS).split(DATA_DELIM);

            for(String sampleMapping : sampleMappings)
            {
                String[] mappingItems = sampleMapping.split(SUB_ITEM_DELIM);
                SourceSampleIds.put(mappingItems[0], mappingItems[1]);
            }
        }

        DriverGenes = Sets.newHashSet();

        if(cmd.hasOption(DRIVER_GENE_PANEL_OPTION))
        {
            try
            {
                DriverGeneFile.read(cmd.getOptionValue(DRIVER_GENE_PANEL_OPTION)).forEach(x -> DriverGenes.add(x.gene()));
            }
            catch(IOException e)
            {
                CMP_LOGGER.error("failed to load driver gene panel file: {}", e.toString());
            }
        }

        if(Categories.containsKey(FUSION) || Categories.containsKey(DRIVER))
        {
            GeneMapping = new GeneNameMapping();
        }
        else
        {
            GeneMapping = null;
        }

    }

    public String sourceSampleId(final String source, final String sampleId)
    {
        String suffix = SourceSampleIds.get(source);
        return suffix != null ? sampleId + suffix : sampleId;
    }

    public boolean isValid() { return mIsValid; }
    public boolean singleSample() { return SampleIds.size() == 1; }
    public boolean multiSample() { return SampleIds.size() > 1; }

    public String getGeneMappedName(final String geneName)
    {
        if(GeneMapping == null || GeneMapping.hasNewGene(geneName))
            return geneName;

        return GeneMapping.getNewName(geneName);
    }

    private void loadSampleIds(final CommandLine cmd)
    {
        if(cmd.hasOption(SAMPLE_ID_FILE))
        {
            SampleIds.addAll(ConfigUtils.loadSampleIdsFile(cmd.getOptionValue(SAMPLE_ID_FILE)));

            if(SampleIds.isEmpty())
            {
                mIsValid = false;
            }
        }
        else if(cmd.hasOption(SAMPLE))
        {
            SampleIds.add(cmd.getOptionValue(SAMPLE));
        }
    }

    private void loadDatabaseSources(final CommandLine cmd)
    {
        if(!cmd.hasOption(DB_SOURCES))
            return;

        // form DB1;db_url;db_user;db_pass DB2;db_url;db_user;db_pass etc
        String dbSourcesStr = cmd.getOptionValue(DB_SOURCES, "");
        String[] dbSources = dbSourcesStr.split(DATA_DELIM, -1);

        if(dbSources.length != 2)
        {
            CMP_LOGGER.error("invalid DB source config({}) - must contain 2 entries", dbSourcesStr);
            mIsValid = false;
            return;
        }

        for(String dbSourceStr : dbSources)
        {
            String[] dbItems = dbSourceStr.split(ITEM_DELIM, -1);

            if(dbItems.length != 4)
            {
                CMP_LOGGER.error("invalid DB source config({})", dbSourceStr);
                mIsValid = false;
                return;
            }

            String sourceName = dbItems[0];
            String dbUrl = "jdbc:" + dbItems[1] + DB_DEFAULT_ARGS;
            String dbUsername = dbItems[2];
            String dbPass = dbItems[3];

            try
            {
                final DatabaseAccess dbAccess = new DatabaseAccess(dbUsername, dbPass, dbUrl);

                if(DbConnections.containsKey(sourceName))
                {
                    CMP_LOGGER.error("repeated DB source name({})", sourceName);
                    mIsValid = false;
                    return;
                }

                DbConnections.put(sourceName, dbAccess);
                SourceNames.add(sourceName);
            }
            catch(SQLException e)
            {
                mIsValid = false;
            }
        }
    }

    private void loadFileSources(final CommandLine cmd)
    {
        if(!cmd.hasOption(FILE_SOURCES))
            return;

        // form: sample_dir=pipe_v1;/path_to_sample_dir/;linx_dir=linx;purple_dir=purple etc OR
        // form: linx_dir=pipe_v1;/path_to_linx_data/;purple_dir/path_to_purple_data/ etc OR
        String fileSourcesStr = cmd.getOptionValue(FILE_SOURCES, "");
        String[] fileSourceEntries = fileSourcesStr.split(DATA_DELIM, -1);

        if(fileSourceEntries.length != 2)
        {
            CMP_LOGGER.error("invalid file source config({}) - must contain 2 entries", fileSourceEntries);
            mIsValid = false;
            return;
        }

        for(String fileSourceStr : fileSourceEntries)
        {
            FileSources fileSources = fromConfig(fileSourceStr);

            if(fileSources == null)
            {
                FileSources.clear();
                return;
            }

            FileSources.put(fileSources.Source, fileSources);
            SourceNames.add(fileSources.Source);
        }
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(
                CATEGORIES, true,
                "Categories to check separated by ';' from: PURITY, DRIVER, SOMATIC_VARIANT, GERMLINE_VARIANT, FUSION, DISRUPTION");

        options.addOption(MATCH_LEVEL, true, "Match level from REPORTABLE (default) or DETAILED");
        options.addOption(SAMPLE, true, "Sample data file");
        options.addOption(SAMPLE_ID_FILE, true, "Sample data file");
        options.addOption(SOURCE_SAMPLE_MAPPINGS, true, "Optional specific source suffixes");
        options.addOption(DRIVER_GENE_PANEL_OPTION, true, DRIVER_GENE_PANEL_OPTION_DESC);
        options.addOption(THRESHOLDS, true, "In form: Field,AbsoluteDiff,PercentDiff, separated by ';'");

        options.addOption(DB_SOURCES, true, "Database configurations keyed by soure name");
        options.addOption(FILE_SOURCES, true, "File locations keyed by source name");
        options.addOption(WRITE_DETAILED_FILES, false, "Write per-type details files");
        options.addOption(THREADS, true, "Thread count (default 0, not multi-threaded)");

        addDatabaseCmdLineArgs(options);
        addOutputOptions(options);
        addLoggingOptions(options);
    }
}
