package com.hartwig.hmftools.common.utils.config;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;

public class ConfigUtils
{
    public static final String LOG_DEBUG = "log_debug";
    public static final String LOG_DEBUG_DESC = "Log at DEBUG level";

    public static final String LOG_LEVEL = "log_level";
    public static final String LOG_LEVEL_DESC = "Specify log level: ERROR, WARN, INFO, DEBUG or TRACE";

    public static final String SAMPLE_ID_FILE = "sample_id_file";
    public static final String SAMPLE_ID_FILE_DESC = "Sample ID CSV file with 'SampleId' column";

    public static final String GENE_ID_FILE = "gene_id_file";
    public static final String GENE_ID_FILE_DESC = "Restricted set of Gene IDs in CSV file";

    public static final String IGNORE_SAMPLE_ID = "#";

    public static final String CONFIG_FILE_DELIM = ",";

    private static final Logger LOGGER = LogManager.getLogger(ConfigUtils.class);

    public static void addLoggingOptions(final ConfigBuilder configBuilder)
    {
        configBuilder.addFlag(LOG_DEBUG, LOG_DEBUG_DESC);
        configBuilder.addConfigItem(LOG_LEVEL, LOG_LEVEL_DESC);
    }

    public static void setLogLevel(final ConfigBuilder configBuilder)
    {
        if(configBuilder.isRegistered(LOG_DEBUG) && configBuilder.hasFlag(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }
        else if(configBuilder.isRegistered(LOG_LEVEL) && configBuilder.hasValue(LOG_LEVEL))
        {
            Configurator.setRootLevel(Level.valueOf(configBuilder.getValue(LOG_LEVEL)));
        }
    }

    public static void addSampleIdFile(final ConfigBuilder configBuilder, boolean required)
    {
        configBuilder.addPath(SAMPLE_ID_FILE, required, SAMPLE_ID_FILE_DESC);
    }

    public static List<String> loadSampleIdsFile(final ConfigBuilder configBuilder)
    {
        if(!configBuilder.hasValue(SAMPLE_ID_FILE))
            return Lists.newArrayList();

        return loadSampleIdsFile(configBuilder.getValue(SAMPLE_ID_FILE));
    }

    public static List<String> loadSampleIdsFile(final String filename)
    {
        return loadDelimitedIdFile(filename, FLD_SAMPLE_ID, CSV_DELIM);
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
            return Collections.emptyList();
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
            if(line.isEmpty() || line.startsWith(IGNORE_SAMPLE_ID))
                continue;

            String[] items = line.split(delim, -1);
            ids.add(items[idIndex]);
        }

        return ids;
    }

    public static String convertWildcardSamplePath(final String samplePath, final String sampleId)
    {
        if(samplePath == null)
            return samplePath;

        return samplePath.replaceAll("\\*", sampleId);
    }

}
