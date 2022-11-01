package com.hartwig.hmftools.compar;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION_DESC;
import static com.hartwig.hmftools.common.utils.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.compar.Category.ALL_CATEGORIES;
import static com.hartwig.hmftools.compar.Category.LINX_CATEGORIES;
import static com.hartwig.hmftools.compar.Category.PURPLE_CATEGORIES;
import static com.hartwig.hmftools.compar.Category.purpleCategories;
import static com.hartwig.hmftools.compar.Category.linxCategories;
import static com.hartwig.hmftools.compar.CommonUtils.DATA_DELIM;
import static com.hartwig.hmftools.compar.CommonUtils.ITEM_DELIM;
import static com.hartwig.hmftools.compar.FileSources.fromConfig;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.DB_DEFAULT_ARGS;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
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

    public final Set<String> DriverGenes;
    public final Set<String> AlternateTranscriptDriverGenes;

    public final DiffThresholds Thresholds;

    public final String OutputDir;
    public final String OutputId;

    public final boolean WriteDetailed;
    public final int Threads;

    private final Map<String,SampleIdMapping> mSampleIdMappings; // if required, mapping from original to new sampleId
    private boolean mIsValid;

    // config strings
    public static final String CATEGORIES = "categories";
    public static final String MATCH_LEVEL = "match_level";

    public static final String DB_SOURCE = "db_source";
    public static final String FILE_SOURCE = "file_source";
    public static final String THRESHOLDS = "thresholds";

    public static final String SAMPLE = "sample";
    public static final String WRITE_DETAILED_FILES = "write_detailed";

    public static final Logger CMP_LOGGER = LogManager.getLogger(ComparConfig.class);

    public static final String REF_SOURCE = "ref";
    public static final String NEW_SOURCE = "new";

    public ComparConfig(final CommandLine cmd)
    {
        mIsValid = true;

        SampleIds = Lists.newArrayList();
        mSampleIdMappings = Maps.newHashMap();

        Categories = Maps.newHashMap();

        MatchLevel matchLevel = MatchLevel.valueOf(cmd.getOptionValue(MATCH_LEVEL, REPORTABLE.toString()));

        String categoriesStr = cmd.getOptionValue(CATEGORIES, ALL_CATEGORIES);

        CMP_LOGGER.info("default match level {}, categories: {}", matchLevel, categoriesStr);

        if(categoriesStr.equals(ALL_CATEGORIES))
        {
            Arrays.stream(Category.values()).forEach(x -> Categories.put(x, matchLevel));
        }
        else if(categoriesStr.contains(PURPLE_CATEGORIES) || categoriesStr.contains(LINX_CATEGORIES))
        {
            if(categoriesStr.contains(PURPLE_CATEGORIES))
                purpleCategories().forEach(x -> Categories.put(x, matchLevel));

            if(categoriesStr.contains(LINX_CATEGORIES))
                linxCategories(). forEach(x -> Categories.put(x, matchLevel));
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

                    CMP_LOGGER.info("specific category({}) and matchLevel({})", category, specificMatchLevel);
                }
                else
                {
                    category = Category.valueOf(catData);
                    specificMatchLevel = matchLevel;
                }

                Categories.put(category, specificMatchLevel);
            }
        }

        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID);
        WriteDetailed = cmd.hasOption(WRITE_DETAILED_FILES);
        Threads = parseThreads(cmd);

        SourceNames = Lists.newArrayList(REF_SOURCE, NEW_SOURCE);
        loadSampleIds(cmd);

        DbConnections = Maps.newHashMap();
        FileSources = Maps.newHashMap();

        if(cmd.hasOption(formConfigSourceStr(DB_SOURCE, REF_SOURCE)) && cmd.hasOption(formConfigSourceStr(DB_SOURCE, NEW_SOURCE)))
        {
            loadDatabaseSources(cmd);
        }
        else if(cmd.hasOption(formConfigSourceStr(FILE_SOURCE, REF_SOURCE)) && cmd.hasOption(formConfigSourceStr(FILE_SOURCE, NEW_SOURCE)))
        {
            loadFileSources(cmd);
        }
        else
        {
            CMP_LOGGER.error("missing DB or file source ref and new config");
            mIsValid = false;
        }

        Thresholds = new DiffThresholds();
        Thresholds.loadConfig(cmd.getOptionValue(THRESHOLDS, ""));

        DriverGenes = Sets.newHashSet();
        AlternateTranscriptDriverGenes = Sets.newHashSet();

        if(cmd.hasOption(DRIVER_GENE_PANEL_OPTION))
        {
            try
            {
                List<DriverGene> driverGenes = DriverGeneFile.read(cmd.getOptionValue(DRIVER_GENE_PANEL_OPTION));

                for(DriverGene driverGene : driverGenes)
                {
                    DriverGenes.add(driverGene.gene());

                    if(!driverGene.additionalReportedTranscripts().isEmpty())
                        AlternateTranscriptDriverGenes.add(driverGene.gene());
                }
            }
            catch(IOException e)
            {
                CMP_LOGGER.error("failed to load driver gene panel file: {}", e.toString());
            }
        }
    }

    public String sourceSampleId(final String source, final String sampleId)
    {
        SampleIdMapping mapping = mSampleIdMappings.get(sampleId);

        if(mapping == null || !mapping.SourceMapping.containsKey(source))
            return sampleId;

        return mapping.SourceMapping.get(source);
    }

    public boolean isValid() { return mIsValid; }
    public boolean singleSample() { return SampleIds.size() == 1; }
    public boolean multiSample() { return SampleIds.size() > 1; }

    private class SampleIdMapping
    {
        public final String SampleId;
        public Map<String,String> SourceMapping;

        public SampleIdMapping(final String sampleId)
        {
            SampleId = sampleId;
            SourceMapping = Maps.newHashMap();
        }
    }

    private static final String COL_SAMPLE_ID = "SampleId";
    private static final String COL_REF_SAMPLE_ID = "RefSampleId";
    private static final String COL_NEW_SAMPLE_ID = "NewSampleId";

    private void loadSampleIds(final CommandLine cmd)
    {
        if(cmd.hasOption(SAMPLE))
        {
            SampleIds.add(cmd.getOptionValue(SAMPLE));
            return;
        }

        if(!cmd.hasOption(SAMPLE_ID_FILE))
        {
            CMP_LOGGER.error("missing sample_id_file or sample config");
            mIsValid = false;
            return;
        }

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(cmd.getOptionValue(SAMPLE_ID_FILE)));
            String header = lines.get(0);
            lines.remove(0);

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DATA_DELIM);

            int sampleIndex = fieldsIndexMap.get(COL_SAMPLE_ID);
            Integer refSampleIndex = fieldsIndexMap.get(COL_REF_SAMPLE_ID);
            Integer newSampleIndex = fieldsIndexMap.get(COL_NEW_SAMPLE_ID);

            for(String line : lines)
            {
                String[] values = line.split(DATA_DELIM, -1);

                String sampleId = values[sampleIndex];

                if(sampleId.startsWith("#"))
                    continue;

                SampleIds.add(sampleId);

                String refSampleId = refSampleIndex != null ? values[refSampleIndex] : null;
                String newSampleId = newSampleIndex != null ? values[newSampleIndex] : null;

                if(refSampleId != null || newSampleId != null)
                {
                    SampleIdMapping mapping = new SampleIdMapping(sampleId);
                    mSampleIdMappings.put(sampleId, mapping);

                    if(refSampleId != null && SourceNames.size() >= 1);
                        mapping.SourceMapping.put(SourceNames.get(0), refSampleId);

                    if(newSampleId != null && SourceNames.size() >= 2)
                        mapping.SourceMapping.put(SourceNames.get(1), newSampleId);
                }
            }

            CMP_LOGGER.info("loaded {} samples from file", SampleIds.size());
        }
        catch(IOException e)
        {
            CMP_LOGGER.error("failed to load sample IDs: {}", e.toString());
        }
    }

    private static String formConfigSourceStr(final String sourceType, final String sourceName)
    {
        return format("%s_%s", sourceType, sourceName);
    }

    private void loadDatabaseSources(final CommandLine cmd)
    {
        if(!cmd.hasOption(formConfigSourceStr(DB_SOURCE, REF_SOURCE)) || !cmd.hasOption(formConfigSourceStr(DB_SOURCE, NEW_SOURCE)))
            return;

        // form DB1;db_url;db_user;db_pass DB2;db_url;db_user;db_pass etc

        for(String sourceName : SourceNames)
        {
            String dbConfigValue = cmd.getOptionValue(formConfigSourceStr(DB_SOURCE, sourceName));
            String[] dbItems = dbConfigValue.split(ITEM_DELIM, -1);

            if(dbItems.length != 3)
            {
                CMP_LOGGER.error("invalid DB source config({})", dbConfigValue);
                mIsValid = false;
                return;
            }

            String dbUrl = "jdbc:" + dbItems[0] + DB_DEFAULT_ARGS;
            String dbUsername = dbItems[1];
            String dbPass = dbItems[2];

            try
            {
                final DatabaseAccess dbAccess = new DatabaseAccess(dbUsername, dbPass, dbUrl);
                DbConnections.put(sourceName, dbAccess);
            }
            catch(SQLException e)
            {
                mIsValid = false;
            }
        }
    }

    private void loadFileSources(final CommandLine cmd)
    {
        // form: sample_dir=pipe_v1;/path_to_sample_dir/;linx_dir=linx;purple_dir=purple etc OR
        // form: linx_dir=pipe_v1;/path_to_linx_data/;purple_dir/path_to_purple_data/ etc OR

        for(String sourceName : SourceNames)
        {
            String fileConfigValue = cmd.getOptionValue(formConfigSourceStr(FILE_SOURCE, sourceName));

            FileSources fileSources = fromConfig(sourceName, fileConfigValue);

            if(fileSources == null)
            {
                FileSources.clear();
                return;
            }

            FileSources.put(fileSources.Source, fileSources);
        }
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(
                CATEGORIES, true,
                "Categories to check separated by ';' from: PURITY, DRIVER, SOMATIC_VARIANT, GERMLINE_VARIANT, FUSION, DISRUPTION");

        options.addOption(MATCH_LEVEL, true, "Match level from REPORTABLE (default) or DETAILED");
        options.addOption(SAMPLE, true, "Sample data file");
        addSampleIdFile(options);
        options.addOption(DRIVER_GENE_PANEL_OPTION, true, DRIVER_GENE_PANEL_OPTION_DESC);
        options.addOption(THRESHOLDS, true, "In form: Field,AbsoluteDiff,PercentDiff, separated by ';'");

        options.addOption(formConfigSourceStr(DB_SOURCE, REF_SOURCE), true, "Database configurations for reference data");
        options.addOption(formConfigSourceStr(DB_SOURCE, NEW_SOURCE), true, "Database configurations for new data");
        options.addOption(formConfigSourceStr(FILE_SOURCE, REF_SOURCE), true, "File locations for reference data");
        options.addOption(formConfigSourceStr(FILE_SOURCE, NEW_SOURCE), true, "File locations for new data");
        options.addOption(WRITE_DETAILED_FILES, false, "Write per-type details files");
        addThreadOptions(options);

        addDatabaseCmdLineArgs(options);
        addOutputOptions(options);
        addLoggingOptions(options);
    }

    public ComparConfig()
    {
        mIsValid = true;

        SampleIds = Lists.newArrayList();
        Categories = Maps.newHashMap();
        OutputDir = null;
        OutputId = "";
        WriteDetailed = false;
        Threads = 0;

        DbConnections = Maps.newHashMap();
        FileSources = Maps.newHashMap();
        SourceNames = Lists.newArrayList();

        Thresholds = new DiffThresholds();
        DriverGenes = Sets.newHashSet();
        AlternateTranscriptDriverGenes = Sets.newHashSet();
        mSampleIdMappings = Maps.newHashMap();
    }
}
