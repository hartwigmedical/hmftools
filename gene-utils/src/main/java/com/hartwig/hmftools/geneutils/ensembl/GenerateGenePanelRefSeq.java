package com.hartwig.hmftools.geneutils.ensembl;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.genepanel.HmfTranscriptRegionFile.DEFAULT_DELIM;
import static com.hartwig.hmftools.common.genome.genepanel.HmfTranscriptRegionFile.EXON_DATA_DELIM;
import static com.hartwig.hmftools.common.genome.genepanel.HmfTranscriptRegionFile.ITEM_DELIM;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.readQueryString;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.writeRecordsAsTsv;
import static com.hartwig.hmftools.geneutils.ensembl.EnsemblDAO.createEnsemblDbConnection;

import java.io.BufferedWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.StringJoiner;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.genepanel.HmfTranscriptRegionFile;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jooq.DSLContext;
import org.jooq.Record;
import org.jooq.Result;

import org.jooq.types.UInteger;
import org.jooq.types.ULong;

public class GenerateGenePanelRefSeq
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

        String geneFile = String.format("%sall_genes.%s.tsv", outputDir, refGenomeVersion.identifier());
        generateGenes(context, geneFile, refGenomeVersion, "sql/ensembl_gene_panel.sql");

        String reqSeqFile = String.format("%sref_seq.%s.tsv", outputDir, refGenomeVersion.identifier());
        generateRefSeqMapping(context, reqSeqFile);

        GU_LOGGER.info("complete");
    }

    private static void generateGenes(
            final DSLContext context, final String outputFile, final RefGenomeVersion refGenomeVersion, String genePanelQuery)
    {
        try
        {
            final Result<Record> results = context.fetch(readQueryString(Resources.getResource(genePanelQuery)));
            GU_LOGGER.debug("gene query returned {} entries", results.size());

            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write(HmfTranscriptRegionFile.header());
            // writer.write(",EntrezIdsHgnc");
            writer.newLine();

            int recordIndex = 0;

            while(recordIndex < results.size())
            {
                Record record = results.get(recordIndex);

                String geneId = (String) record.get("gene_id");

                String chromosome = (String) record.get("chromosome");
                String geneName = (String) record.get("gene_name");
                UInteger geneStart = (UInteger) record.get("gene_start");
                UInteger geneEnd = (UInteger) record.get("gene_end");
                Byte strand = (Byte) record.get("strand");

                String entrezIds = record.get("entrezId") == null ? "" : (String) record.get("entrezId");
                // String entrezIdsHgnc = record.get("entrezId") == null ? "" : (String) record.get("entrezIdHgnc");
                String chrBand = (String) record.get("chromosome_band");

                String transId = (String) record.get("transcript_id");
                UInteger transStart = (UInteger) record.get("transcript_start");
                UInteger transEnd = (UInteger) record.get("transcript_end");
                ULong codingStart = (ULong) record.get("coding_start");
                ULong codingEnd = (ULong) record.get("coding_end");
                StringJoiner exonData = new StringJoiner(ITEM_DELIM);

                while(recordIndex < results.size())
                {
                    Record nextRecord = results.get(recordIndex);
                    String nextTransId = (String) nextRecord.get("transcript_id");

                    if(!transId.equals(nextTransId))
                        break;

                    UInteger exonStart = (UInteger) nextRecord.get("exon_start");
                    UInteger exonEnd = (UInteger) nextRecord.get("exon_end");
                    exonData.add(String.format("%d%s%d", exonStart.intValue(), EXON_DATA_DELIM, exonEnd.intValue()));
                    ++recordIndex;
                }

                StringJoiner sj = new StringJoiner(DEFAULT_DELIM);
                sj.add(geneId);
                sj.add(geneName);
                sj.add(refGenomeVersion.versionedChromosome(chromosome));
                sj.add(String.valueOf(geneStart));
                sj.add(String.valueOf(geneEnd));
                sj.add(String.valueOf(strand));
                sj.add(entrezIds.replaceAll(",", ITEM_DELIM));
                sj.add(chrBand);

                sj.add(transId);
                sj.add(String.valueOf(transStart));
                sj.add(String.valueOf(transEnd));
                sj.add(codingStart != null ? String.valueOf(codingStart) : "");
                sj.add(codingEnd != null ? String.valueOf(codingEnd) : "");
                sj.add(exonData.toString());

                writer.write(sj.toString());
                // writer.write(String.format(",%s", entrezIdsHgnc));
                writer.newLine();
            }

            writer.close();

            GU_LOGGER.info("written gene panel output to {}", outputFile);
        }
        catch (final IOException e)
        {
            GU_LOGGER.error("error writing Ensembl gene data file: {}", e.toString());
            return;
        }
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
