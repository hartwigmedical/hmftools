package com.hartwig.hmftools.genepanel;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.google.common.io.Resources;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
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

    private static final String OUT_PATH = "out";
    private static final String ENSEMBL_VERSION = "ensembl";

    private static final String ENSEMBLDB_URL_37 = "jdbc:mysql://ensembldb.ensembl.org:3337/homo_sapiens_core_89_37";
    private static final String ENSEMBLDB_URL_38 = "jdbc:mysql://ensembldb.ensembl.org:3306/homo_sapiens_core_89_38";
    private static final String DB_USER = "anonymous";

    public static void main(String[] args) throws ParseException, IOException, SQLException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String outputFilePath = cmd.getOptionValue(OUT_PATH);
        final String ensemblVersion = cmd.getOptionValue(ENSEMBL_VERSION);

        if (outputFilePath == null || !VERSIONS.contains(ensemblVersion)) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("HmfGenePanelBuilder", options);
            System.exit(1);
        }

        final String database = ensemblVersion.equals("37") ? ENSEMBLDB_URL_37 : ENSEMBLDB_URL_38;

        LOGGER.info("Querying " + database);
        final Result<Record> queryResults = queryEnsembldb(database);
        writeFile(cmd, queryResults);
        LOGGER.info("Written output to " + new File(outputFilePath).getAbsolutePath());
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(OUT_PATH, true, "Path towards the tsv output file.");
        options.addOption(ENSEMBL_VERSION, true, "Ensembl version to use. Must be either 37 or 38.");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static String readEnsemblQuery() throws IOException {
        final List<String> lines = Resources.readLines(Resources.getResource("ensembl_query.sql"), Charset.defaultCharset());
        return StringUtils.join(lines.toArray(), "\n");
    }

    @NotNull
    private static Result<Record> queryEnsembldb(@NotNull final String database) throws SQLException, IOException {
        // MIVO: disable annoying jooq self-ad message
        System.setProperty("org.jooq.no-logo", "true");
        final Connection conn = DriverManager.getConnection(database, DB_USER, "");
        final DSLContext context = DSL.using(conn, SQLDialect.MYSQL);
        final String query = readEnsemblQuery();
        return context.fetch(query);
    }

    private static void writeFile(@NotNull final CommandLine cmd, @NotNull final Result<Record> records) throws IOException {
        final BufferedWriter writer = new BufferedWriter(new FileWriter(cmd.getOptionValue(OUT_PATH), false));
        // MIVO: format as tsv without header containing column names
        final CSVFormat format = new CSVFormat().header(false).delimiter('\t').nullString("").quoteString("");
        writer.write(records.formatCSV(format));
        writer.close();
    }
}
