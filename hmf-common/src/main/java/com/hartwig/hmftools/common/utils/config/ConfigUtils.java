package com.hartwig.hmftools.common.utils.config;

import static com.hartwig.hmftools.common.utils.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class ConfigUtils
{
    public static final String LOG_DEBUG = "log_debug";
    public static final String LOG_DEBUG_DESC = "Log at DEBUG level";

    public static final String LOG_LEVEL = "log_level";
    public static final String LOG_LEVEL_DESC = "Specify log level: ERROR, WARN, INFO, DEBUG or TRACE";

    public static final String SAMPLE_ID_FILE = "sample_id_file";
    public static final String SAMPLE_ID_FILE_DESC = "Sample ID CSV file with 'SampleId' column";

    public static final String SAMPLE_ID_COLUMN = "SampleId";

    private static final Logger LOGGER = LogManager.getLogger(ConfigUtils.class);

    public static void addLoggingOptions(final Options options)
    {
        options.addOption(LOG_DEBUG, false, LOG_DEBUG_DESC);
        options.addOption(LOG_LEVEL, true, LOG_LEVEL_DESC);
    }

    public static void addLoggingOptions(final ConfigBuilder configBuilder)
    {
        configBuilder.addFlagItem(LOG_DEBUG, LOG_DEBUG_DESC);
        configBuilder.addConfigItem(LOG_LEVEL, LOG_LEVEL_DESC);
    }

    public static void setLogLevel(final CommandLine cmd)
    {
        if(cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }
        else if(cmd.hasOption(LOG_LEVEL))
        {
            Configurator.setRootLevel(Level.valueOf(cmd.getOptionValue(LOG_LEVEL)));
        }
    }

    public static void setLogLevel(final ConfigBuilder configBuilder)
    {
        if(configBuilder.hasFlag(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }
        else if(configBuilder.hasValue(LOG_LEVEL))
        {
            Configurator.setRootLevel(Level.valueOf(configBuilder.getValue(LOG_LEVEL)));
        }
    }

    public static void addSampleIdFile(final Options options)
    {
        options.addOption(SAMPLE_ID_FILE, true, SAMPLE_ID_FILE_DESC);
    }

    public static void addSampleIdFile(final ConfigBuilder configBuilder, boolean required)
    {
        configBuilder.addPathItem(SAMPLE_ID_FILE, required, SAMPLE_ID_FILE_DESC);
    }

    public static double getConfigValue(@NotNull final CommandLine cmd, final String configName, double defaultValue)
    {
        return cmd.hasOption(configName) ? Double.parseDouble(cmd.getOptionValue(configName)) : defaultValue;
    }

    public static int getConfigValue(@NotNull final CommandLine cmd, final String configName, int defaultValue)
    {
        return cmd.hasOption(configName) ? Integer.parseInt(cmd.getOptionValue(configName)) : defaultValue;
    }

    public static boolean getConfigValue(@NotNull final CommandLine cmd, final String configName, boolean defaultValue)
    {
        return cmd.hasOption(configName) ? Boolean.parseBoolean(cmd.getOptionValue(configName)) : defaultValue;
    }

    public static List<String> loadSampleIdsFile(final ConfigBuilder configBuilder)
    {
        if(!configBuilder.hasValue(SAMPLE_ID_FILE))
            return Lists.newArrayList();

        return loadSampleIdsFile(configBuilder.getValue(SAMPLE_ID_FILE));
    }

    public static List<String> loadSampleIdsFile(final CommandLine cmd)
    {
        if(!cmd.hasOption(SAMPLE_ID_FILE))
            return Lists.newArrayList();

        return loadSampleIdsFile(cmd.getOptionValue(SAMPLE_ID_FILE));
    }

    public static List<String> loadSampleIdsFile(final String filename)
    {
        return loadDelimitedIdFile(filename, SAMPLE_ID_COLUMN, CSV_DELIM);
    }

    public static List<String> loadGeneIdsFile(final String filename)
    {
        return loadDelimitedIdFile(filename, "GeneId", CSV_DELIM);
    }

    public static List<String> loadDelimitedIdFile(final String filename, final String idColumn, final String delim)
    {
        // the file must either contain the idColumn in the header, or have a single column
        final List<String> ids = Lists.newArrayList();

        if(filename == null || filename.isEmpty())
        {
            LOGGER.error("{} file({}) missing or invalid", idColumn, filename);
            return ids;
        }

        try
        {
            List<String> fileContents = Files.readAllLines(new File(filename).toPath());
            return loadDelimitedIdFile(fileContents, idColumn, delim);
        }
        catch (IOException e)
        {
            LOGGER.error("failed to read {} file({}): {}", idColumn, filename, e.toString());
            return Collections.EMPTY_LIST;
        }
    }

    @VisibleForTesting
    public static List<String> loadDelimitedIdFile(final List<String> fileContents, final String idColumn, final String delim)
    {
        final List<String> ids = Lists.newArrayList();

        String firstLine = fileContents.get(0);
        int columnCount = firstLine.split(delim, -1).length;
        int idIndex = -1;

        if(firstLine.contains(idColumn))
        {
            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(firstLine, delim);
            idIndex = fieldsIndexMap.get(idColumn);
            fileContents.remove(0);
        }
        else
        {
            if(columnCount > 1)
                return ids;

            idIndex = 0;
        }

        for(String line : fileContents)
        {
            if(line.startsWith("#") || line.isEmpty())
                continue;

            String[] items = line.split(delim, -1);
            ids.add(items[idIndex]);
        }

        return ids;
    }

    public static boolean containsFlag(final CommandLine cmd, final String opt)
    {
        if(cmd.hasOption(opt))
        {
            LOGGER.info("Using non default {} flag", opt);
            return true;
        }
        return false;
    }

}
