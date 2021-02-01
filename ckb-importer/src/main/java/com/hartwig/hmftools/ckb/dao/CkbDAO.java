package com.hartwig.hmftools.ckb.dao;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;

import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrial;
import com.hartwig.hmftools.ckb.datamodel.drug.Drug;

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

public class CkbDAO {
    private static final Logger LOGGER = LogManager.getLogger(CkbDAO.class);

    private static final String DEV_CATALOG = "ckb_test";

    @NotNull
    private final DSLContext context;

    @NotNull
    public static CkbDAO connectToCkbDAO(@NotNull final String userName, @NotNull final String password, @NotNull final String url)
            throws SQLException {
        final Connection conn = DriverManager.getConnection(url, userName, password);
        final String catalog = conn.getCatalog();
        LOGGER.info("Connecting to database {}", catalog);

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
      //  ClinicalTrialDAO.clearClinicalTrial(context);
       // DrugDAO.clearDrug(context);
    }

    public void writeCkb(@NotNull CkbEntry ckbEntry) {
        int count = 0;
        LOGGER.info("Starting writing clinical trial");
        for (ClinicalTrial clinicalTrial : ckbEntry.clinicalTrial()) {
            //ClinicalTrialDAO.writeClinicalTrial(context, clinicalTrial);
            count = counting(count, "clinical trial object", ckbEntry.clinicalTrial().size());
        }
        LOGGER.info("Finished writing clinical trial");

        LOGGER.info("Starting writing drug object");
        count = 0;
        for (Drug drug : ckbEntry.drug()) {
           // DrugDAO.writeDrug(context, drug);
            count = counting(count, "drug object", ckbEntry.drug().size());
        }
        LOGGER.info("Finished writing drug object");

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
