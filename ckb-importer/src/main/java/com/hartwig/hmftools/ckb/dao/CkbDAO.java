package com.hartwig.hmftools.ckb.dao;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;

import com.sun.scenario.Settings;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CkbDAO {
    private static final Logger LOGGER = LogManager.getLogger(CkbDAO.class);

    private static final String DEV_CATALOG = "ckb_test";



//    @NotNull
//    private final DSLContext context;
//
//    @NotNull
//    public static CkbDAO connectToCkbDAO(@NotNull final String userName, @NotNull final String password, @NotNull final String url)
//            throws SQLException {
//        final Connection conn = DriverManager.getConnection(url, userName, password);
//        final String catalog = conn.getCatalog();
//        LOGGER.info("Connecting to database {}", catalog);
//
//        return new CkbDAO(DSL.using(conn, SQLDialect.MYSQL, settings(catalog)));
//    }
//
//    @Nullable
//    private static Settings settings(@NotNull String catalog) {
//        if (catalog.equals(DEV_CATALOG)) {
//            return null;
//        }
//
//        return new Settings().withRenderMapping(new RenderMapping().withSchemata(new MappedSchema().withInput(DEV_CATALOG)
//                .withOutput(catalog)));
//    }
//
//    private CkbDAO(@NotNull final DSLContext context) {
//        this.context = context;
//    }







}
