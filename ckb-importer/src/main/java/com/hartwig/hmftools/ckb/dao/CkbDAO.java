package com.hartwig.hmftools.ckb.dao;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;

import com.hartwig.hmftools.ckb.datamodel.CkbJsonDatabase;

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

    public void deleteAll() {

    }

    public void writeCkb(@NotNull CkbJsonDatabase ckbEntry) {

    }

    private int counting(int count, @NotNull String specificObject, int totalEntriesOfObject) {
        count++;
        if (count % 1000 == 0) {
            LOGGER.info(" Completed inserting {} of {} CKB entries into CKB db of the {} entries",
                    count,
                    specificObject,
                    totalEntriesOfObject);
        }
        return count;
    }
}
