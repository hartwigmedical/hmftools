package com.hartwig.hmftools.serve.dao;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.common.serve.actionability.ActionableEvents;
import com.hartwig.hmftools.serve.extraction.KnownEvents;
import com.hartwig.hmftools.serve.extraction.events.EventInterpretation;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.DSLContext;
import org.jooq.SQLDialect;
import org.jooq.conf.MappedSchema;
import org.jooq.conf.RenderMapping;
import org.jooq.conf.Settings;
import org.jooq.impl.DSL;

public class ServeDatabaseAccess {

    private static final Logger LOGGER = LogManager.getLogger(ServeDatabaseAccess.class);
    private static final String DEV_CATALOG = "serve_test";

    public static final String DB_USER = "db_user";
    public static final String DB_PASS = "db_pass";
    public static final String DB_URL = "db_url";

    public static final String DB_DEFAULT_ARGS = "?serverTimezone=UTC&useSSL=false";

    @NotNull
    private final Connection connection;
    @NotNull
    private final DSLContext context;

    @NotNull
    private final ServeDAO serveDAO;

    public ServeDatabaseAccess(@NotNull final String userName, @NotNull final String password, @NotNull final String url)
            throws SQLException {
        System.setProperty("org.jooq.no-logo", "true");
        System.setProperty("org.jooq.no-tips", "true");

        this.connection = DriverManager.getConnection(url, userName, password);
        String catalog = connection.getCatalog();
        LOGGER.debug("Connecting to database '{}'", catalog);
        this.context = DSL.using(connection, SQLDialect.MYSQL, settings(catalog));

        this.serveDAO = new ServeDAO(context);
    }

    @Nullable
    private static Settings settings(@NotNull String catalog) {
        if (catalog.equals(DEV_CATALOG)) {
            return null;
        }

        return new Settings().withRenderMapping(new RenderMapping().withSchemata(new MappedSchema().withInput(DEV_CATALOG)
                .withOutput(catalog)));
    }

    @NotNull
    public static ServeDatabaseAccess databaseAccess(@NotNull CommandLine cmd) throws SQLException {
        return databaseAccess(cmd, false);
    }

    @NotNull
    public static ServeDatabaseAccess databaseAccess(@NotNull CommandLine cmd, boolean applyDefaultArgs) throws SQLException {
        String userName = cmd.getOptionValue(DB_USER);
        String password = cmd.getOptionValue(DB_PASS);
        String databaseUrl = cmd.getOptionValue(DB_URL);
        String jdbcUrl = "jdbc:" + databaseUrl;

        if (applyDefaultArgs && !jdbcUrl.contains("serverTimezone") && !jdbcUrl.contains("useSSL")) {
            jdbcUrl += DB_DEFAULT_ARGS;
        }

        return new ServeDatabaseAccess(userName, password, jdbcUrl);
    }

    public static void addDatabaseCmdLineArgs(@NotNull Options options) {
        addDatabaseCmdLineArgs(options, false);
    }

    public static void addDatabaseCmdLineArgs(@NotNull Options options, boolean isRequired) {
        options.addOption(Option.builder(DB_USER).desc("Database username").hasArg(true).required(isRequired).build());
        options.addOption(Option.builder(DB_PASS).desc("Database password").hasArg(true).required(isRequired).build());
        options.addOption(Option.builder(DB_URL).desc("Database url").hasArg(true).required(isRequired).build());
    }

    public void writeServeDAO(@NotNull ActionableEvents actionableEvents, KnownEvents knownEvents,
            List<EventInterpretation> eventInterpretations) {
        serveDAO.write(actionableEvents, knownEvents, eventInterpretations);
    }
}
