package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.DEVTAG;
import static com.hartwig.hmftools.vicc.database.Tables.FEATURENAME;
import static com.hartwig.hmftools.vicc.database.Tables.GENE;
import static com.hartwig.hmftools.vicc.database.Tables.GENEIDENTIFIER;
import static com.hartwig.hmftools.vicc.database.Tables.TAG;
import static com.hartwig.hmftools.vicc.database.Tables.VICCENTRY;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.vicc.ViccJsonToSQLImporter;
import com.hartwig.hmftools.vicc.datamodel.GeneIdentifier;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
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

        writeTags(id, viccEntry.tags());
        writeDevTags(id, viccEntry.devTags());
        writeGeneIdentifiers(id, viccEntry.geneIdentifiers());
        writeGenes(id, viccEntry.genes());
        writeFeatureNames(id, viccEntry.featureNames());
    }

    private void writeTags(int viccEntryId, @NotNull List<String> tags) {
        for (String tag : tags) {
            context.insertInto(TAG, TAG.TAGNAME, TAG.VICCENTRYID).values(tag, viccEntryId).execute();
        }
    }

    private void writeDevTags(int viccEntryId, @NotNull List<String> devTags) {
        for (String devTag : devTags) {
            context.insertInto(DEVTAG, DEVTAG.DEVTAGNAME, DEVTAG.VICCENTRYID).values(devTag, viccEntryId).execute();
        }
    }

    private void writeGeneIdentifiers(int viccEntryId, @NotNull List<GeneIdentifier> geneIdentifiers) {
        for (GeneIdentifier geneIdentifier : geneIdentifiers) {
            context.insertInto(GENEIDENTIFIER,
                    GENEIDENTIFIER.SYMBOL,
                    GENEIDENTIFIER.ENTREZID,
                    GENEIDENTIFIER.ENSEMBLGENEID,
                    GENEIDENTIFIER.VICCENTRYID)
                    .values(geneIdentifier.symbol(), geneIdentifier.entrezId(), geneIdentifier.ensemblGeneId(), viccEntryId)
                    .execute();
        }
    }

    private void writeGenes(int viccEntryId, @NotNull List<String> genes) {
        for (String gene : genes) {
            context.insertInto(GENE, GENE.GENENAME, GENE.VICCENTRYID).values(gene, viccEntryId).execute();
        }
    }

    private void writeFeatureNames(int viccEntryId, @Nullable List<String> featureNames) {
        if (featureNames != null) {
            for (String featureName : featureNames) {
                context.insertInto(FEATURENAME, FEATURENAME.FEATURENAME_, FEATURENAME.VICCENTRYID)
                        .values(featureName, viccEntryId)
                        .execute();
            }
        }
    }

    public void deleteAll() {
        LOGGER.info("Deleting all from vicc db");

        context.deleteFrom(TAG).execute();
        context.deleteFrom(DEVTAG).execute();
        context.deleteFrom(GENEIDENTIFIER).execute();
        context.deleteFrom(GENE).execute();
        context.deleteFrom(FEATURENAME).execute();
        context.deleteFrom(VICCENTRY).execute();
    }
}
