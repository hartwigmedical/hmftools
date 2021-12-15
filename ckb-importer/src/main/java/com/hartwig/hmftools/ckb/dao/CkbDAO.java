package com.hartwig.hmftools.ckb.dao;

import static com.hartwig.hmftools.ckb.database.tables.Ckbentry.CKBENTRY;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;

import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrial;
import com.hartwig.hmftools.ckb.datamodel.evidence.Evidence;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;

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
    private final TherapyDAO therapyDAO;
    @NotNull
    private final IndicationDAO indicationDAO;
    @NotNull
    private final VariantDAO variantDAO;
    @NotNull
    private final EvidenceDAO evidenceDAO;
    @NotNull
    private final ClinicalTrialDAO clinicalTrialDAO;

    @NotNull
    public static CkbDAO connectToCkbDAO(@NotNull String userName, @NotNull String password, @NotNull String url) throws SQLException {
        LOGGER.info("Connecting to database CKB at {}", url);
        Connection conn = DriverManager.getConnection(url, userName, password);
        String catalog = conn.getCatalog();

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
        this.therapyDAO = new TherapyDAO(context);
        this.indicationDAO = new IndicationDAO(context);
        this.variantDAO = new VariantDAO(context);
        this.evidenceDAO = new EvidenceDAO(context, therapyDAO, indicationDAO);
        this.clinicalTrialDAO = new ClinicalTrialDAO(context, therapyDAO, indicationDAO);
    }

    public void deleteAll() {
        // Note that deletions should go from branch to root
        variantDAO.deleteAll();

        // The evidence and trial DAOs need to emptied first before therapy and indication, to avoid key violations.
        evidenceDAO.deleteAll();
        clinicalTrialDAO.deleteAll();

        therapyDAO.deleteAll();
        indicationDAO.deleteAll();

        context.deleteFrom(CKBENTRY).execute();
    }

    public void write(@NotNull CkbEntry ckbEntry) {
        int id = context.insertInto(CKBENTRY, CKBENTRY.CKBPROFILEID, CKBENTRY.CREATEDATE, CKBENTRY.UPDATEDATE, CKBENTRY.PROFILENAME)
                .values(ckbEntry.profileId(),
                        ckbEntry.createDate(),
                        ckbEntry.updateDate(),
                        ckbEntry.profileName())
                .returning(CKBENTRY.ID)
                .fetchOne()
                .getValue(CKBENTRY.ID);

        for (Variant variant : ckbEntry.variants()) {
            variantDAO.write(variant, id);
        }

        for (Evidence evidence : ckbEntry.evidences()) {
            evidenceDAO.write(evidence, id);
        }

        for (ClinicalTrial clinicalTrial : ckbEntry.clinicalTrials()) {
            clinicalTrialDAO.write(clinicalTrial, id);
        }
    }
}
