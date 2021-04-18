package com.hartwig.hmftools.ensembl;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URL;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jooq.CSVFormat;
import org.jooq.DSLContext;
import org.jooq.Record;
import org.jooq.Result;
import org.jooq.SQLDialect;
import org.jooq.impl.DSL;
import org.jooq.tools.StringUtils;

public class HmfEnsemblResourceBuilder {

    private static final Logger LOGGER = LogManager.getLogger(HmfEnsemblResourceBuilder.class);

    private static final String GENE_TSV = "gene_tsv";
    private static final String REFSEQ_MAPPING_TSV = "refseq_mapping_tsv";

    private static final String ENSEMBL_DB_URL_37 = "jdbc:mysql://ensembldb.ensembl.org:3337/homo_sapiens_core_89_37";
    private static final String ENSEMBL_DB_URL_38 = "jdbc:mysql://ensembldb.ensembl.org:3306/homo_sapiens_core_89_38";
    private static final String ENSEMBL_DB_USER = "anonymous";

    public static void main(String[] args) throws ParseException, IOException, SQLException {
        final Options options = createOptions();
        final CommandLine cmd = new DefaultParser().parse(options, args);

        final RefGenomeVersion refGenomeVersion = RefGenomeVersion.from(cmd.getOptionValue(RefGenomeVersion.REF_GENOME_VERSION));
        final String geneTsv = cmd.getOptionValue(GENE_TSV);
        final String refseqMappingTsv = cmd.getOptionValue(REFSEQ_MAPPING_TSV);

        if (geneTsv == null || refseqMappingTsv == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("HmfEnsemblResourceBuilder", options);
            System.exit(1);
        }

        final String database = refGenomeVersion.is37() ? ENSEMBL_DB_URL_37 : ENSEMBL_DB_URL_38;

        LOGGER.info("Running hmf-ensembl-resource-builder with version {}", refGenomeVersion);
        final DSLContext context = connectEnsemblDb(database);
        LOGGER.info(" Connected to {}", database);

        generateGenes(context, geneTsv, refGenomeVersion);
        generateRefSeqMapping(context, refseqMappingTsv);

        LOGGER.info("Complete");
    }

    private static void generateGenes(@NotNull final DSLContext context, @NotNull final String outputTsv,
            @NotNull final RefGenomeVersion refGenomeVersion) throws IOException {
        final Result<Record> geneResult = context.fetch(read(Resources.getResource("sql/ensembl_gene_query.sql")));
        LOGGER.info(" Gene query returned {} entries", geneResult.size());

        writeAsTsv(outputTsv, geneResult);
        LOGGER.info(" Written gene output to {}", outputTsv);

        if (refGenomeVersion == RefGenomeVersion.V38) {
            LOGGER.info("  Patching genes to v38 in {}", outputTsv);
            final List<String> genes = Files.readAllLines(new File(outputTsv).toPath());
            final List<String> updatedGenes = Lists.newArrayList();

            // First column is the chromosome, and we need to prefix with "chr" to patch the genes to v38 ref genome
            for (String gene : genes) {
                updatedGenes.add("chr" + gene);
            }
            Files.write(new File(outputTsv).toPath(), updatedGenes);
            LOGGER.info("  Finished patching genes to v38 in {}", outputTsv);
        }
    }

    private static void generateRefSeqMapping(@NotNull final DSLContext context, @NotNull final String outputTsv) throws IOException {
        final Result<Record> refseqMappingResult = context.fetch(read(Resources.getResource("sql/ensembl_refseq_mapping_query.sql")));
        LOGGER.info(" RefSeq mapping query returned {} entries", refseqMappingResult.size());

        writeAsTsv(outputTsv, refseqMappingResult);
        LOGGER.info(" Written RefSeq mapping output to {}", outputTsv);
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(RefGenomeVersion.REF_GENOME_VERSION, true, "Ref genome version to generate files for");
        options.addOption(GENE_TSV, true, "Path towards the gene tsv output file.");
        options.addOption(REFSEQ_MAPPING_TSV, true, "Path towards the refseq mapping tsv output file.");
        return options;
    }

    @NotNull
    private static DSLContext connectEnsemblDb(@NotNull final String database) throws SQLException {
        // Disable annoying jooq self-ad message
        System.setProperty("org.jooq.no-logo", "true");
        final Connection conn = DriverManager.getConnection(database, ENSEMBL_DB_USER, "");
        return DSL.using(conn, SQLDialect.MYSQL);
    }

    @NotNull
    private static String read(@NotNull final URL queryResource) throws IOException {
        final List<String> lines = Resources.readLines(queryResource, Charset.defaultCharset());
        return StringUtils.join(lines.toArray(), "\n");
    }

    private static void writeAsTsv(@NotNull final String tsv, @NotNull final Result<Record> records) throws IOException {
        final BufferedWriter writer = new BufferedWriter(new FileWriter(tsv, false));
        // Format as tsv without header containing column names
        final CSVFormat format = new CSVFormat().header(false).delimiter('\t').nullString("").quoteString("");
        writer.write(records.formatCSV(format));
        writer.close();
    }
}
