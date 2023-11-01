package com.hartwig.hmftools.geneutils.ensembl;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.readQueryString;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.writeRecordsAsTsv;
import static com.hartwig.hmftools.geneutils.ensembl.EnsemblDAO.createEnsemblDbConnection;

import java.io.IOException;
import java.sql.SQLException;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

import org.jooq.DSLContext;
import org.jooq.Record;
import org.jooq.Result;

public class RunAdHocQuery
{
    private static final String QUERY_FILE = "query_file";
    private static final String OUTPUT_FILE = "output_file";

    public static void main(String[] args) throws IOException, SQLException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addConfigItem(REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC, V37.toString());
        configBuilder.addPath(QUERY_FILE, true, "Ad-hoc query SQL file");
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output file for test query results");

        EnsemblDAO.addCmdLineArgs(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);
        addOutputDir(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        RefGenomeVersion refGenomeVersion = RefGenomeVersion.from(configBuilder);
        String outputFile = configBuilder.getValue(OUTPUT_FILE);

        if(refGenomeVersion == null || outputFile == null)
        {
            GU_LOGGER.error("missing config");
            System.exit(1);
        }

        DSLContext context = createEnsemblDbConnection(configBuilder);

        if(context == null)
            System.exit(1);

        GU_LOGGER.debug("database connection established");

        String query = readQueryString(configBuilder.getValue(QUERY_FILE));

        final Result<Record> results = context.fetch(query);

        GU_LOGGER.info("query returned {} entries", results.size());
        writeRecordsAsTsv(outputFile, results);

        GU_LOGGER.info("complete");
    }
}
