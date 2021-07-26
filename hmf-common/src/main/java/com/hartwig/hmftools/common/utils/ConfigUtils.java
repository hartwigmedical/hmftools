package com.hartwig.hmftools.common.utils;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;

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

    private static final Logger LOGGER = LogManager.getLogger(ConfigUtils.class);

    public static void addLoggingOptions(final Options options)
    {
        options.addOption(LOG_DEBUG, false, "Log verbose");
        options.addOption(LOG_LEVEL, true, "Specify log level");
    }

    public static void setLogLevel(final CommandLine cmd)
    {
        if (cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }
        else if(cmd.hasOption(LOG_LEVEL))
        {
            Configurator.setRootLevel(Level.valueOf(cmd.getOptionValue(LOG_LEVEL)));
        }
    }

    public static double getConfigValue(final CommandLine cmd, final String configName, double defaultValue)
    {
        return cmd.hasOption(configName) ? Double.parseDouble(cmd.getOptionValue(configName)) : defaultValue;
    }

    public static int getConfigValue(final CommandLine cmd, final String configName, int defaultValue)
    {
        return cmd.hasOption(configName) ? Integer.parseInt(cmd.getOptionValue(configName)) : defaultValue;
    }

    public static boolean getConfigValue(final CommandLine cmd, final String configName, boolean defaultValue)
    {
        return cmd.hasOption(configName) ? Boolean.parseBoolean(cmd.getOptionValue(configName)) : defaultValue;
    }

    public static List<String> loadSampleIdFile(final String filename)
    {
        final List<String> sampleIds = Lists.newArrayList();

        if(filename == null || filename.isEmpty())
        {
            LOGGER.error("sampleId file({}) missing or invalid", filename);
            return sampleIds;
        }

        try
        {
            Files.readAllLines(new File(filename).toPath()).stream()
                    .filter(x -> !x.equalsIgnoreCase("SampleId"))
                    .filter(x -> !x.isEmpty())
                    .forEach(x -> sampleIds.add(x));
        }
        catch (IOException e)
        {
            LOGGER.error("failed to read sampleId file({}): {}", filename, e.toString());
        }

        return sampleIds;
    }

    public static List<String> loadDelimitedSampleIdFile(final String filename, final String delim)
    {
        final List<String> sampleIds = Lists.newArrayList();

        if(filename == null || filename.isEmpty())
        {
            LOGGER.error("sampleId file({}) missing or invalid", filename);
            return sampleIds;
        }

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            if(fileContents.isEmpty())
                return sampleIds;

            String header = fileContents.get(0);

            if(!header.contains("SampleId"))
                return sampleIds;

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, delim);
            fileContents.remove(0);
            int sampleIdIndex = fieldsIndexMap.get("SampleId");

            for(String sampleData : fileContents)
            {
                String[] items = sampleData.split(delim, -1);
                sampleIds.add(items[sampleIdIndex]);
            }
        }
        catch (IOException e)
        {
            LOGGER.error("failed to read sampleId file({}): {}", filename, e.toString());
        }

        return sampleIds;
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

    @NotNull
    public static String readableFile(@NotNull final CommandLine cmd, @NotNull final String opt) throws IOException, ParseException
    {
        if(!cmd.hasOption(opt))
        {
            throw new ParseException(opt + " is a required option");
        }
        final String file = cmd.getOptionValue(opt);
        if(!new File(file).exists())
        {
            throw new IOException("Unable to read file:  " + file);
        }

        return file;
    }

    @NotNull
    public static String writableOutputDirectory(@NotNull final CommandLine cmd, @NotNull final String opt)
            throws ParseException, IOException
    {
        if(!cmd.hasOption(opt))
        {
            throw new ParseException(opt + " is a required option");
        }

        final String outputDirString = cmd.getOptionValue(opt);
        final File outputDir = new File(outputDirString);
        if(!outputDir.exists() && !outputDir.mkdirs())
        {
            throw new IOException("Unable to write output directory " + outputDirString);
        }

        return outputDirString;
    }

}
