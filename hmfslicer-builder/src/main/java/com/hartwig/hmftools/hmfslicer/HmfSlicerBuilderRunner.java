package com.hartwig.hmftools.hmfslicer;

import java.io.File;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.io.reader.FileReader;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
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

public final class HmfSlicerBuilderRunner {
    private static final Logger LOGGER = LogManager.getLogger(HmfSlicerBuilderRunner.class);

    private static final String FASTQ_FILE = "file";
    private static final String FASTQ_ROOT_DIR = "dir";
    private static final String CSV_OUT_DIR = "out";
    private static final String THREAD_COUNT = "threadCount";
    private static final String ENSEMBLDB_URL = "ensembldb.ensembl.org";
    private static final String DB_USER = "anonymous";

    public static void main(String[] args)
            throws ParseException, IOException, InterruptedException, SQLException, EmptyFileException {
        final Connection conn = DriverManager.getConnection(ENSEMBLDB_URL, DB_USER, "");
        final DSLContext context = DSL.using(conn, SQLDialect.MYSQL);
        final String baseQuery = readEnsemblQuery();
        final String genes = readGeneList();
        final String groupByClause = "group by gene_name";
        final String finalQuery = baseQuery + "and display_label in (" + genes + ") " + groupByClause;
        final Result<Record> queryResults = context.fetch(finalQuery);
        for (final Record queryResult : queryResults) {
            System.out.println(queryResult.get(0));
        }
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(FASTQ_FILE, true, "Path towards the original fastq file.");
        options.addOption(FASTQ_ROOT_DIR, true, "Path towards the flowcell dir.");
        options.addOption(CSV_OUT_DIR, true, "Path towards the csv output file.");
        options.addOption(THREAD_COUNT, true, "Number of max threads to use (only used when running on a directory).");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    private static String readGeneList() throws IOException, EmptyFileException {
        final List<String> lines = FileReader.build().readLines(
                new File(Resources.getResource("gene_panel").getPath()).toPath());
        return StringUtils.join(lines, ",");
    }

    private static String readEnsemblQuery() throws IOException, EmptyFileException {
        final List<String> lines = FileReader.build().readLines(
                new File(Resources.getResource("ensembl_query.sql").getPath()).toPath());
        return StringUtils.join(lines, "\n");
    }
}
