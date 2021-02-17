package com.hartwig.hmftools.ckb.dao;

import static com.hartwig.hmftools.ckb.database.tables.Ckbentry.CKBENTRY;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.time.LocalDate;

import com.hartwig.hmftools.ckb.datamodel.CkbEntry;

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

public final class CkbDAO {

    private static final Logger LOGGER = LogManager.getLogger(CkbDAO.class);

    private static final String DEV_CATALOG = "ckb_test";

    @NotNull
    private final DSLContext context;

    @NotNull
    public static CkbDAO connectToCkbDAO(@NotNull String userName, @NotNull String password, @NotNull String url) throws SQLException {
        Connection conn = DriverManager.getConnection(url, userName, password);
        String catalog = conn.getCatalog();
        LOGGER.info("Connecting to database CKB {}", catalog);

        return new CkbDAO(DSL.using(conn, SQLDialect.MYSQL, settings(catalog)));
    }

    @Nullable
    private static org.jooq.conf.Settings settings(@NotNull String catalog) {
        if (catalog.equals(DEV_CATALOG)) {
            return null;
        }

        return new Settings().withRenderMapping(new RenderMapping().withSchemata(new MappedSchema().withInput(DEV_CATALOG)
                .withOutput(catalog)));
    }

    private CkbDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    public void write(@NotNull CkbEntry ckbEntry) {
        int id = context.insertInto(CKBENTRY, CKBENTRY.CKBPROFILEID,
                CKBENTRY.PROFILENAME,
                CKBENTRY.CREATEDATE,
                CKBENTRY.UPDATEDATE)
                .values(ckbEntry.profileId(),
                        ckbEntry.profileName(),
                        sqlDate(ckbEntry.createDate()),
                        sqlDate(ckbEntry.updateDate()))
                .returning(CKBENTRY.ID)
                .fetchOne()
                .getValue(CKBENTRY.ID);
    }

    public void deleteAll() {

    }

    @Nullable
    private static java.sql.Date sqlDate(@Nullable LocalDate date) {
        return date != null ? java.sql.Date.valueOf(date) : null;
    }
}
