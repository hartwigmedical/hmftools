package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.VICCENTRY;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;

import com.hartwig.hmftools.vicc.ViccJsonToSQLImporter;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.SQLDialect;
import org.jooq.impl.DSL;

public class ViccDAO {

    private static final Logger LOGGER = LogManager.getLogger(ViccJsonToSQLImporter.class);

    @NotNull
    private final DSLContext context;

    public static ViccDAO connectToViccDAO(@NotNull final String userName, @NotNull final String password, @NotNull final String url)
            throws SQLException {
        final Connection conn = DriverManager.getConnection(url, userName, password);
        final String catalog = conn.getCatalog();

        LOGGER.debug("Connecting to database {}", catalog);
        return new ViccDAO(DSL.using(conn, SQLDialect.MYSQL));
    }

    private ViccDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    public void writeViccEntry(@NotNull ViccEntry viccEntry) {
        int id = context.insertInto(VICCENTRY, VICCENTRY.SOURCE)
                .values(viccEntry.source())
                .returning(VICCENTRY.ID)
                .fetchOne()
                .getValue(VICCENTRY.ID);

        LOGGER.info("Added VICC entry with id " + id + ": " + viccEntry);
    }

    public void deleteAll() {
        LOGGER.info("Deleting all from vicc db");
        context.deleteFrom(VICCENTRY).execute();
    }
}
