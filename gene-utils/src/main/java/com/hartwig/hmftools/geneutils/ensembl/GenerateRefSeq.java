package com.hartwig.hmftools.geneutils.ensembl;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.readQueryString;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.writeRecordsAsTsv;
import static com.hartwig.hmftools.geneutils.ensembl.EnsemblDAO.createEnsemblDbConnection;

import java.io.IOException;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.Map;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

import org.apache.commons.cli.ParseException;
import org.jooq.DSLContext;
import org.jooq.Record;
import org.jooq.Result;

public class GenerateRefSeq
{
    private static final String OUTPUT_FILE = "output_file";

    public static void main(String[] args) throws ParseException, IOException, SQLException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addConfigItem(REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC);
        EnsemblDAO.addCmdLineArgs(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);

        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output file for RefSeq transcript mappings");

        configBuilder.checkAndParseCommandLine(args);

        RefGenomeVersion refGenomeVersion = RefGenomeVersion.from(configBuilder);

        GU_LOGGER.info("writing Ensembl RefSeq transcript mapping, ref-genome-version({})", refGenomeVersion);

        DSLContext context = createEnsemblDbConnection(configBuilder);

        if(context == null)
            System.exit(1);

        GU_LOGGER.debug("database connection established");

        String reqSeqOutputFile = configBuilder.getValue(OUTPUT_FILE);
        generateRefSeqMapping(context, reqSeqOutputFile);

        GU_LOGGER.info("ReqSeq transcript mapping complete");
    }

    static Map<String, String> getRefSeqMapping(final DSLContext context)
    {
        Map<String, String> result = new HashMap<>();
        final Result<Record> refseqMappingResult = context.fetch(readQueryString(Resources.getResource("ensembl_sql/ensembl_refseq_mapping.sql")));
        GU_LOGGER.debug("RefSeq mapping query returned {} entries", refseqMappingResult.size());
        refseqMappingResult.forEach(record -> {
            String refSeqId = (String) record.get("refSeqId");
            String transcriptId = (String) record.get("transcriptId");
            result.put(transcriptId, refSeqId);
        });
        return result;
    }

    private static void generateRefSeqMapping(final DSLContext context, final String outputFile)
    {
        final Result<Record> refseqMappingResult = context.fetch(readQueryString(Resources.getResource("ensembl_sql/ensembl_refseq_mapping.sql")));
        GU_LOGGER.debug("RefSeq mapping query returned {} entries", refseqMappingResult.size());

        writeRecordsAsTsv(outputFile, refseqMappingResult);
        GU_LOGGER.info("written RefSeq mapping output to {}", outputFile);
    }
}
