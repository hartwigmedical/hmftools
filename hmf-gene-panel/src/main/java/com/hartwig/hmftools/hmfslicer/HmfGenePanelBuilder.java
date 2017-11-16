package com.hartwig.hmftools.hmfslicer;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.EmptyFileException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.Record;
import org.jooq.Result;
import org.jooq.SQLDialect;
import org.jooq.impl.DSL;
import org.jooq.tools.StringUtils;

public class HmfGenePanelBuilder {

    private static final Logger LOGGER = LogManager.getLogger(HmfGenePanelBuilder.class);

    private static final String OUT_PATH = "out";
    private static final String DATABASE = "homo_sapiens_core_89_37";
    private static final String ENSEMBLDB_URL = "jdbc:mysql://ensembldb.ensembl.org:3337/" + DATABASE;
    private static final String DB_USER = "anonymous";

    public static void main(String[] args) throws ParseException, IOException, InterruptedException, SQLException, EmptyFileException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String outputFilePath = cmd.getOptionValue(OUT_PATH);
        if (outputFilePath == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("HmfGenePanelBuilder", options);
        } else {
            final Result<Record> queryResults = queryEnsembldb();
            writeFile(cmd, queryResults);
            LOGGER.info("Written output to " + new File(outputFilePath).getAbsolutePath());
        }
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(OUT_PATH, true, "Path towards the tsv output file.");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static String readEnsemblQuery() throws IOException, EmptyFileException {
        final List<String> lines = Resources.readLines(Resources.getResource("ensembl_query.sql"), Charset.defaultCharset());
        return StringUtils.join(lines.toArray(), "\n");
    }

    @NotNull
    @VisibleForTesting
    static Result<Record> queryEnsembldb() throws SQLException, IOException, EmptyFileException {
        // MIVO: disable annoying jooq self-ad message
        System.setProperty("org.jooq.no-logo", "true");
        final Connection conn = DriverManager.getConnection(ENSEMBLDB_URL, DB_USER, "");
        final DSLContext context = DSL.using(conn, SQLDialect.MYSQL);
        final String query = readEnsemblQuery();
        return context.fetch(query);
    }

    private static void writeFile(@NotNull final CommandLine cmd, @NotNull final Result<Record> records) throws IOException {
        final BufferedWriter writer = new BufferedWriter(new FileWriter(cmd.getOptionValue(OUT_PATH), false));
        // MIVO: format as tsv without header containing column names
        writer.write(records.formatCSV(false, '\t', ""));
        writer.close();
    }
}
