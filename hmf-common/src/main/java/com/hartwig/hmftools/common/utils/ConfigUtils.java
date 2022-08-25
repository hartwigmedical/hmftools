package com.hartwig.hmftools.common.utils;

import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class ConfigUtils
{
    public static final String LOG_DEBUG = "log_debug";
    public static final String LOG_LEVEL = "log_level";
    public static final String CSV_DELIM = ",";

    private static final Logger LOGGER = LogManager.getLogger(ConfigUtils.class);

    public static void addLoggingOptions(final Options options)
    {
        options.addOption(LOG_DEBUG, false, "Log at DEBUG level");
        options.addOption(LOG_LEVEL, true, "Specify log level: ERROR, WARN, INFO, DEBUG or TRACE");
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

    public static List<String> loadSampleIdsFile(final String filename)
    {
        return loadDelimitedIdFile(filename, "SampleId", CSV_DELIM);
    }

    public static List<String> loadGeneIdsFile(final String filename)
    {
        return loadDelimitedIdFile(filename, "GeneId", CSV_DELIM);
    }

    public static List<String> loadDelimitedIdFile(final String filename, final String idColumn, final String delim)
    {
        final List<String> ids = Lists.newArrayList();

        if(filename == null || filename.isEmpty())
        {
            LOGGER.error("{} file({}) missing or invalid", idColumn, filename);
            return ids;
        }

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            if(fileContents.isEmpty())
                return ids;

            String header = fileContents.get(0);

            if(!header.contains(idColumn))
                return ids;

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, delim);
            fileContents.remove(0);
            int idIndex = fieldsIndexMap.get(idColumn);

            for(String line : fileContents)
            {
                if(line.startsWith("#") || line.isEmpty())
                    continue;

                String[] items = line.split(delim, -1);
                ids.add(items[idIndex]);
            }
        }
        catch (IOException e)
        {
            LOGGER.error("failed to read {} file({}): {}", idColumn, filename, e.toString());
        }

        return ids;
    }

    @NotNull
    public static <E extends Enum<E>> E defaultEnumValue(@NotNull final CommandLine cmd, @NotNull final String argument,
            @NotNull final E defaultValue) throws ParseException
    {
        if(cmd.hasOption(argument))
        {
            final String optionValue = cmd.getOptionValue(argument);
            try
            {
                final E value = E.valueOf(defaultValue.getDeclaringClass(), optionValue);
                if(!value.equals(defaultValue))
                {
                    LOGGER.info("Using non default value {} for parameter {}", optionValue, argument);
                }

                return value;
            } catch(IllegalArgumentException e)
            {
                throw new ParseException("Invalid validation stringency: " + optionValue);
            }
        }

        return defaultValue;
    }

    public static boolean containsFlag(@NotNull final CommandLine cmd, @NotNull final String opt)
    {
        if(cmd.hasOption(opt))
        {
            LOGGER.info("Using non default {} flag", opt);
            return true;
        }
        return false;
    }

}
