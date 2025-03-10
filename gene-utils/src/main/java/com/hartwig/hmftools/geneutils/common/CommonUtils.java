package com.hartwig.hmftools.geneutils.common;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jooq.CSVFormat;
import org.jooq.Record;
import org.jooq.Result;
import org.jooq.tools.StringUtils;

public final class CommonUtils
{
    public static final String APP_NAME = "GeneUtils";

    public static final Logger GU_LOGGER = LogManager.getLogger(CommonUtils.class);

    public static final String RESOURCE_REPO_DIR = "resource_repo_dir";
    public static final String RESOURCE_REPO_DIR_DESC = "The directory holding the public HMF resources repo";
    public static final String ENSEMBL_DIR = "ensembl_data_cache";
    public static final String DRIVER_GENE_PANEL_TSV = "driver_gene_panel";

    public static String getEnsemblDirectory(final RefGenomeVersion refGenomeVersion, final String resourceRepoDir)
    {
        return resourceRepoDir + ENSEMBL_DIR + File.separator + refGenomeVersion.identifier();
    }

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

            final CSVFormat format = new CSVFormat().header(true).delimiter(TSV_DELIM).nullString("").quoteString("");
            writer.write(records.formatCSV(format));
            writer.close();
        }
        catch (final IOException e)
        {
            GU_LOGGER.error("error writing query results file: {}", e.toString());
        }
    }

    public static boolean createOutputDir(final String outputDir)
    {
        final File dir = new File(outputDir);
        if(!dir.exists() && !dir.mkdirs())
        {
            GU_LOGGER.error("unable to write directory " + outputDir);
            return false;
        }

        return true;
    }
}
