package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.PMKB;
import static com.hartwig.hmftools.vicc.database.Tables.PMKBGENE;
import static com.hartwig.hmftools.vicc.database.Tables.PMKBTISSUE;
import static com.hartwig.hmftools.vicc.database.Tables.PMKBTUMOR;
import static com.hartwig.hmftools.vicc.database.Tables.PMKBVARIANT;

import com.hartwig.hmftools.vicc.datamodel.pmkb.Pmkb;
import com.hartwig.hmftools.vicc.datamodel.pmkb.PmkbTissue;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

final class PmkbDAOFunctions {

    private PmkbDAOFunctions() {
    }

    static void write(@NotNull DSLContext context, int viccEntryId, @NotNull Pmkb pmkb) {
        int id = context.insertInto(PMKB, PMKB.VICCENTRYID).values(viccEntryId).returning(PMKB.ID).fetchOne().getValue(PMKB.ID);

        context.insertInto(PMKBTUMOR, PMKBTUMOR.NAME, PMKBTUMOR.IDTUMOR, PMKBTUMOR.PMKBID)
                .values(pmkb.tumor().name(), pmkb.tumor().id(), id)
                .execute();

        for (PmkbTissue tissue : pmkb.tissues()) {
            context.insertInto(PMKBTISSUE, PMKBTISSUE.NAME, PMKBTISSUE.IDTISSUE, PMKBTISSUE.PMKBID)
                    .values(tissue.name(), tissue.id(), id)
                    .execute();
        }

        int variantId = context.insertInto(PMKBVARIANT,
                PMKBVARIANT.NAME,
                PMKBVARIANT.COORDINATES,
                PMKBVARIANT.CHROMOSOME,
                PMKBVARIANT.CYTOBAND,
                PMKBVARIANT.TRANSCRIPT,
                PMKBVARIANT.EFFECT,
                PMKBVARIANT.CODONS,
                PMKBVARIANT.EXONS,
                PMKBVARIANT.DNACHANGE,
                PMKBVARIANT.AMINOACIDCHANGE,
                PMKBVARIANT.GERMLINE,
                PMKBVARIANT.PARTNERGENE,
                PMKBVARIANT.CNVTYPE,
                PMKBVARIANT.CHROMOSOMEBASEDCNV,
                PMKBVARIANT.VARIANTTYPE,
                PMKBVARIANT.COSMIC,
                PMKBVARIANT.DESCRIPTION,
                PMKBVARIANT.DESCRIPTIONTYPE,
                PMKBVARIANT.NOTES,
                PMKBVARIANT.IDVARIANT,
                PMKBVARIANT.PMKBID)
                .values(pmkb.variant().name(),
                        pmkb.variant().coordinates(),
                        pmkb.variant().chromosome(),
                        pmkb.variant().cytoband(),
                        pmkb.variant().transcript(),
                        pmkb.variant().effect(),
                        pmkb.variant().codons(),
                        pmkb.variant().exons(),
                        pmkb.variant().dnaChange(),
                        pmkb.variant().aminoAcidChange(),
                        pmkb.variant().germline(),
                        pmkb.variant().partnerGene(),
                        pmkb.variant().cnvType(),
                        pmkb.variant().chromosomeBasedCnv(),
                        pmkb.variant().variantType(),
                        pmkb.variant().cosmic(),
                        pmkb.variant().description(),
                        pmkb.variant().descriptionType(),
                        pmkb.variant().notes(),
                        pmkb.variant().id(),
                        id)
                .returning(PMKBVARIANT.ID)
                .fetchOne()
                .getValue(PMKBVARIANT.ID);

        context.insertInto(PMKBGENE,
                PMKBGENE.NAME,
                PMKBGENE.CREATEDAT,
                PMKBGENE.UPDATEDAT,
                PMKBGENE.ACTIVEIND,
                PMKBGENE.EXTERNALID,
                PMKBGENE.DESCRIPTION,
                PMKBGENE.IDGENE,
                PMKBGENE.PMKBVARIANTID)
                .values(pmkb.variant().gene().name(),
                        pmkb.variant().gene().createdAt(),
                        pmkb.variant().gene().updatedAt(),
                        pmkb.variant().gene().activeInd(),
                        pmkb.variant().gene().externalId(),
                        pmkb.variant().gene().description(),
                        pmkb.variant().gene().id(),
                        variantId)
                .execute();
    }

    static void deleteAll(@NotNull DSLContext context) {
        context.deleteFrom(PMKBGENE).execute();
        context.deleteFrom(PMKBVARIANT).execute();
        context.deleteFrom(PMKBTISSUE).execute();
        context.deleteFrom(PMKBTUMOR).execute();

        context.deleteFrom(PMKB).execute();
    }
}
