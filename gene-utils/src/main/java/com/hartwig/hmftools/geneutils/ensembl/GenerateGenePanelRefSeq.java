package com.hartwig.hmftools.geneutils.ensembl;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;
import static com.hartwig.hmftools.geneutils.ensembl.EnsemblDAO.createEnsemblDbConnection;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.sql.SQLException;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jooq.DSLContext;
import org.jooq.Record;
import org.jooq.Result;

import org.jooq.CSVFormat;
import org.jooq.tools.StringUtils;
import org.jooq.types.UInteger;

public class GenerateGenePanelRefSeq
{
    private static final String ENSEMBL_DB_URL_37 = "jdbc:mysql://ensembldb.ensembl.org:3337/homo_sapiens_core_89_37";
    private static final String ENSEMBL_DB_URL_38 = "jdbc:mysql://ensembldb.ensembl.org:3306/homo_sapiens_core_89_38";
    private static final String ENSEMBL_DB_USER = "anonymous";

    public static void main(String[] args) throws ParseException, IOException, SQLException
    {
        final Options options = createOptions();
        final CommandLine cmd = new DefaultParser().parse(options, args);

        setLogLevel(cmd);

        RefGenomeVersion refGenomeVersion = RefGenomeVersion.from(cmd.getOptionValue(RefGenomeVersion.REF_GENOME_VERSION));
        String outputDir = parseOutputDir(cmd);
        String geneTsv = String.format("%sall_genes.%s.tsv", outputDir, refGenomeVersion.identifier());
        String refseqMappingTsv = String.format("%sref_seq.%s.tsv", outputDir, refGenomeVersion.identifier());

        if(geneTsv == null || refseqMappingTsv == null)
        {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("HmfEnsemblResourceBuilder", options);
            System.exit(1);
        }

        GU_LOGGER.info("writing Ensembl gene panel and ref-seq data files, ref-genome-version({})", refGenomeVersion);

        DSLContext context = createEnsemblDbConnection(cmd);

        generateGenes(context, geneTsv, refGenomeVersion);
        generateRefSeqMapping(context, refseqMappingTsv);

        GU_LOGGER.info("complete");
    }

    private static void generateGenes(
            final DSLContext context, final String outputTsv, final RefGenomeVersion refGenomeVersion) throws IOException
    {
        final Result<Record> results = context.fetch(read(Resources.getResource("sql/ensembl_gene_query.sql")));
        GU_LOGGER.debug("gene query returned {} entries", results.size());

        writeAsTsv(outputTsv, results);
        GU_LOGGER.info("written gene panel output to {}", outputTsv);

        /*
            if(refGenomeVersion.is38())
            {
                final List<String> genes = Files.readAllLines(new File(outputTsv).toPath());
                final List<String> updatedGenes = Lists.newArrayList();

                // First column is the chromosome, and we need to prefix with "chr" to patch the genes to v38 ref genome
                for(String gene : genes)
                {
                    updatedGenes.add("chr" + gene);
                }
                Files.write(new File(outputTsv).toPath(), updatedGenes);
            }
         */

        // Chromosome	GeneStart	GeneEnd	GeneId	GeneName	EntrezId	ChrBand
        // TranscriptId	TranscriptVersion	TranscriptStart	TranscriptEnd
        // ExonId	ExonStart	ExonEnd	Strand	CodingStart	CodingEnd

        try
        {
            BufferedWriter writer = createBufferedWriter(outputTsv, false);

            StringJoiner header = new StringJoiner("\t");
            header.add("Chromosome").add("GeneStart").add("GeneEnd").add("GeneId").add("GeneName");
            header.add("EntrezId").add("ChrBand");
            header.add("TranscriptId").add("TransVersion").add("TransStart").add("TransEnd");
            header.add("ExonId").add("ExonStart").add("ExonEnd");
            header.add("Strand").add("CodingStart").add("CodingEnd");
            writer.write(header.toString());
            writer.newLine();

            for(final Record record : results)
            {
                StringJoiner sj = new StringJoiner("\t");
                sj.add((String)record.get("chromosome"));
                sj.add(String.valueOf(record.get("gene_start")));
                sj.add(String.valueOf(record.get("gene_end")));
                sj.add((String)record.get("gene_id"));
                sj.add((String)record.get("gene_name"));

                Object entrezId = record.get("entrezId");
                sj.add(entrezId != null ? entrezId.toString() : "");
                sj.add((String)record.get("chromosome_band"));

                sj.add((String)record.get("transcript_id"));
                sj.add(String.valueOf(record.get("transcript_version")));
                sj.add(String.valueOf(record.get("transcript_start")));
                sj.add(String.valueOf(record.get("transcript_end")));

                sj.add((String)record.get("exon_id"));
                sj.add(String.valueOf(record.get("exon_start")));
                sj.add(String.valueOf(record.get("exon_end")));

                sj.add(String.valueOf(record.get("strand")));

                Object codingStart = record.get("coding_start");
                sj.add(codingStart != null ? codingStart.toString() : "");

                Object codingEnd = record.get("coding_end");
                sj.add(codingEnd != null ? codingEnd.toString() : "");

                writer.write(sj.toString());
                writer.newLine();
            }

            writer.close();
        }
        catch (final IOException e)
        {
            GU_LOGGER.error("error writing Ensembl gene data file: {}", e.toString());
        }

    }

    private static void generateRefSeqMapping(final DSLContext context, final String outputTsv) throws IOException
    {
        final Result<Record> refseqMappingResult = context.fetch(read(Resources.getResource("sql/ensembl_refseq_mapping_query.sql")));
        GU_LOGGER.info("RefSeq mapping query returned {} entries", refseqMappingResult.size());

        writeAsTsv(outputTsv, refseqMappingResult);
        GU_LOGGER.info("written RefSeq mapping output to {}", outputTsv);
    }

    private static String read(final URL queryResource) throws IOException
    {
        final List<String> lines = Resources.readLines(queryResource, Charset.defaultCharset());
        return StringUtils.join(lines.toArray(), "\n");
    }

    private static void writeAsTsv(final String outputFile, final Result<Record> records)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            final CSVFormat format = new CSVFormat().header(false).delimiter('\t').nullString("").quoteString("");
            writer.write(records.formatCSV(format));
        }
        catch (final IOException e)
        {
            GU_LOGGER.error("error writing Ensembl gene data file: {}", e.toString());
        }
    }

    private static Options createOptions()
    {
        final Options options = new Options();
        options.addOption(RefGenomeVersion.REF_GENOME_VERSION, true, "Ref genome version (V37 or V38))");
        options.addOption(ENSEMBL_DATA_DIR, true, "Path to Ensembl data cache files");
        EnsemblDAO.addCmdLineArgs(options);
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(LOG_DEBUG, false, "Log verbose");
        return options;
    }

}
