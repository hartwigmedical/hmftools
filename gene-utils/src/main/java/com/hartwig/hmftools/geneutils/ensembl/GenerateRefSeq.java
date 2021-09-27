package com.hartwig.hmftools.geneutils.ensembl;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.readQueryString;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.writeRecordsAsTsv;
import static com.hartwig.hmftools.geneutils.ensembl.EnsemblDAO.createEnsemblDbConnection;

import java.io.IOException;
import java.sql.SQLException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jooq.DSLContext;
import org.jooq.Record;
import org.jooq.Result;

public class GenerateRefSeq
{
    public static void main(String[] args) throws ParseException, IOException, SQLException
    {
        final Options options = createOptions();
        final CommandLine cmd = new DefaultParser().parse(options, args);

        setLogLevel(cmd);

        RefGenomeVersion refGenomeVersion = RefGenomeVersion.from(cmd.getOptionValue(RefGenomeVersion.REF_GENOME_VERSION));
        String outputDir = parseOutputDir(cmd);

        if(refGenomeVersion == null || outputDir == null)
        {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("HmfEnsemblResourceBuilder", options);
            System.exit(1);
        }

        GU_LOGGER.info("writing Ensembl gene panel and ref-seq data files, ref-genome-version({})", refGenomeVersion);

        DSLContext context = createEnsemblDbConnection(cmd);

        if(context == null)
            System.exit(1);

        GU_LOGGER.debug("database connection established");

        String reqSeqFile = String.format("%sref_seq.%s.tsv", outputDir, refGenomeVersion.identifier());
        generateRefSeqMapping(context, reqSeqFile);

        GU_LOGGER.info("complete");
    }

    private static void generateRefSeqMapping(final DSLContext context, final String outputFile)
    {
        final Result<Record> refseqMappingResult = context.fetch(readQueryString(Resources.getResource("sql/ensembl_refseq_mapping.sql")));
        GU_LOGGER.info("RefSeq mapping query returned {} entries", refseqMappingResult.size());

        writeRecordsAsTsv(outputFile, refseqMappingResult);
        GU_LOGGER.info("written RefSeq mapping output to {}", outputFile);
    }

    private static Options createOptions()
    {
        final Options options = new Options();
        options.addOption(RefGenomeVersion.REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC);
        EnsemblDAO.addCmdLineArgs(options);
        addEnsemblDir(options);
        addLoggingOptions(options);
        addOutputDir(options);
        return options;
    }

}
