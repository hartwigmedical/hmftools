package com.hartwig.hmftools.geneutils.common;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.net.URL;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.io.Resources;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jooq.CSVFormat;
import org.jooq.Record;
import org.jooq.Result;
import org.jooq.tools.StringUtils;

public class CommonUtils
{
    public static final Logger GU_LOGGER = LogManager.getLogger(CommonUtils.class);

    public static String readQueryString(final URL queryResource)
    {
        try
        {
            List<String> lines = Resources.readLines(queryResource, Charset.defaultCharset());
            return StringUtils.join(lines.toArray(), "\n");
        }
        catch (final IOException e)
        {
            GU_LOGGER.error("failed to load query resource({}): {}", queryResource, e.toString());
            return "";
        }
    }

    public static String readQueryString(final String filename)
    {
        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            return StringUtils.join(lines.toArray(), "\n");
        }
        catch (final IOException e)
        {
            GU_LOGGER.error("failed to load query file({}): {}", filename, e.toString());
            return "";
        }
    }

    public static void writeRecordsAsTsv(final String outputFile, final Result<Record> records)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            final CSVFormat format = new CSVFormat().header(true).delimiter('\t').nullString("").quoteString("");
            writer.write(records.formatCSV(format));
            writer.close();
        }
        catch (final IOException e)
        {
            GU_LOGGER.error("error writing query results file: {}", e.toString());
        }
    }

}
