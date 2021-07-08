package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.compar.Category.ALL_CATEGORIES;
import static com.hartwig.hmftools.compar.CommonUtils.DATA_DELIM;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.DB_DEFAULT_ARGS;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.ConfigUtils;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ComparConfig
{
    public final List<String> SampleIds;

    public final Map<Category,MatchLevel> Categories;

    // database access
    public final Map<String,DatabaseAccess> DbConnections;
    public final List<String> DbSourceNames;

    public final String OutputDir;

    private boolean mIsValid;

    // config strings
    public static final String CATEGORIES = "categories";
    public static final String MATCH_LEVEL = "match_level";

    public static final String DB_SOURCES = "db_sources";

    public static final String SAMPLE = "sample";
    public static final String SAMPLE_ID_FILE = "sample_id_file";

    public static final Logger CMP_LOGGER = LogManager.getLogger(ComparConfig.class);

    public ComparConfig(final CommandLine cmd)
    {
        mIsValid = true;

        SampleIds = Lists.newArrayList();
        loadSampleIds(cmd);

        Categories = Maps.newHashMap();

        MatchLevel matchLevel = cmd.hasOption(MATCH_LEVEL) ? MatchLevel.valueOf(cmd.getOptionValue(MATCH_LEVEL)) : REPORTABLE;

        if(!cmd.hasOption(CATEGORIES) || cmd.getOptionValue(CATEGORIES).equals(ALL_CATEGORIES))
        {
            Arrays.stream(Category.values()).forEach(x -> Categories.put(x, matchLevel));
        }
        else
        {
            final String[] catDataList = cmd.getOptionValue(CATEGORIES).split(";");

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

        DbConnections = Maps.newHashMap();
        DbSourceNames = Lists.newArrayList();
        loadDatabaseSources(cmd);
    }

    public boolean isValid() { return mIsValid; }

    private void loadSampleIds(final CommandLine cmd)
    {
        if(cmd.hasOption(SAMPLE_ID_FILE))
        {
            SampleIds.addAll(ConfigUtils.loadSampleIdFile(cmd.getOptionValue(SAMPLE_ID_FILE)));

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
        // form DB1;db_url;db_user;db_pass DB2;db_url;db_user;db_pass etc
        String dbSourcesStr = cmd.getOptionValue(DB_SOURCES, "");
        String[] dbSources = dbSourcesStr.split(DATA_DELIM, -1);

        for(String dbSourceStr : dbSources)
        {
            String[] dbItems = dbSourceStr.split(";", -1);

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
                DbSourceNames.add(sourceName);
            }
            catch(SQLException e)
            {
                mIsValid = false;
            }
        }

    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(
                CATEGORIES, true,
                "Categories to check separated by ';' from: DRIVER, LINX_DATA, FUSION, DISRUPTION");

        options.addOption(MATCH_LEVEL, true, "Match level from REPORTABLE, MODERATE or DETAILED");
        options.addOption(SAMPLE, true, "Sample data file");
        options.addOption(SAMPLE_ID_FILE, true, "Sample data file");

        options.addOption(DB_SOURCES, true, "Directory containing standard sample files from pipeline");

        addDatabaseCmdLineArgs(options);

        options.addOption(OUTPUT_DIR, true, "Path to output files");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");
    }
}
