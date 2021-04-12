package com.hartwig.hmftools.genepanel;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URL;
import java.nio.charset.Charset;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.google.common.io.Resources;

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

public class HmfGenePanelBuilder {

    private static final Logger LOGGER = LogManager.getLogger(HmfGenePanelBuilder.class);
    private static final Set<String> VERSIONS = Sets.newHashSet("37", "38");

    private static final String ENSEMBL_VERSION = "ensembl";
    private static final String GENE_TSV = "gene_tsv";
    private static final String REFSEQ_TO_ENSEMBL_TSV = "refseq_to_ensembl_tsv";

    private static final String ENSEMBLDB_URL_37 = "jdbc:mysql://ensembldb.ensembl.org:3337/homo_sapiens_core_89_37";
    private static final String ENSEMBLDB_URL_38 = "jdbc:mysql://ensembldb.ensembl.org:3306/homo_sapiens_core_89_38";
    private static final String DB_USER = "anonymous";

    public static void main(String[] args) throws ParseException, IOException, SQLException {
        final Options options = createOptions();
        final CommandLine cmd = new DefaultParser().parse(options, args);

        final String ensemblVersion = cmd.getOptionValue(ENSEMBL_VERSION);
        final String geneTsv = cmd.getOptionValue(GENE_TSV);
        final String refseqToEnsemblTsv = cmd.getOptionValue(REFSEQ_TO_ENSEMBL_TSV);

        if (!VERSIONS.contains(ensemblVersion) || geneTsv == null || refseqToEnsemblTsv == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("HmfGenePanelBuilder", options);
            System.exit(1);
        }

        final String database = ensemblVersion.equals("37") ? ENSEMBLDB_URL_37 : ENSEMBLDB_URL_38;

        LOGGER.info("Running hmf-gene-panel-builder with version {}", ensemblVersion);
        final DSLContext context = connectEnsemblDb(database);
        LOGGER.info(" Connected to {}", database);

        final Result<Record> geneResult = context.fetch(read(Resources.getResource("sql/ensembl_gene_query.sql")));
        LOGGER.info(" Gene query returned {} entries", geneResult.size());

        writeAsTsv(geneTsv, geneResult);
        LOGGER.info(" Written gene output to {}", geneTsv);

        final Result<Record> refseqToEnsemblResult = context.fetch(read(Resources.getResource("sql/refseq_to_ensembl_transcripts.sql")));
        LOGGER.info(" RefSeq to Ensembl query returned {} entries", refseqToEnsemblResult.size());

        writeAsTsv(refseqToEnsemblTsv, refseqToEnsemblResult);
        LOGGER.info(" Written RefSeq to Ensembl output to {}", refseqToEnsemblTsv);

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(ENSEMBL_VERSION, true, "Ensembl version to use. Must be either 37 or 38.");
        options.addOption(GENE_TSV, true, "Path towards the gene tsv output file.");
        options.addOption(REFSEQ_TO_ENSEMBL_TSV, true, "Path towards the refseq to ensembl tsv output file.");
        return options;
    }

    @NotNull
    private static DSLContext connectEnsemblDb(@NotNull final String database) throws SQLException {
        // Disable annoying jooq self-ad message
        System.setProperty("org.jooq.no-logo", "true");
        final Connection conn = DriverManager.getConnection(database, DB_USER, "");
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
