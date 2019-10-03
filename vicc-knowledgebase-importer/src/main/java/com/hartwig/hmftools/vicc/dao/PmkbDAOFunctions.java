package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.PMKB;
import static com.hartwig.hmftools.vicc.database.Tables.PMKBGENE;
import static com.hartwig.hmftools.vicc.database.Tables.PMKBTISSUE;
import static com.hartwig.hmftools.vicc.database.Tables.PMKBTUMOR;
import static com.hartwig.hmftools.vicc.database.Tables.PMKBVARIANT;

import com.hartwig.hmftools.vicc.datamodel.pmkb.Pmkb;
import com.hartwig.hmftools.vicc.datamodel.pmkb.PmkbGene;
import com.hartwig.hmftools.vicc.datamodel.pmkb.PmkbTissue;
import com.hartwig.hmftools.vicc.datamodel.pmkb.PmkbTumor;
import com.hartwig.hmftools.vicc.datamodel.pmkb.PmkbVariant;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

final class PmkbDAOFunctions {

    private PmkbDAOFunctions() {
    }

    static void write(@NotNull DSLContext context, int viccEntryId, @NotNull Pmkb pmkb) {
        int id = context.insertInto(PMKB, PMKB.VICCENTRYID).values(viccEntryId).returning(PMKB.ID).fetchOne().getValue(PMKB.ID);

        for (PmkbTumor tumor : pmkb.tumor()) {
            context.insertInto(PMKBTUMOR, PMKBTUMOR.IDTUMOR, PMKBTUMOR.NAME, PMKBTUMOR.VICCENTRYID)
                    .values(tumor.id(), tumor.name(), id)
                    .execute();
        }

        for (PmkbTissue tissue : pmkb.tissue()) {
            context.insertInto(PMKBTISSUE, PMKBTISSUE.IDTISSUE, PMKBTISSUE.NAME, PMKBTISSUE.VICCENTRYID)
                    .values(tissue.id(), tissue.name(), id)
                    .execute();
        }

        for (PmkbVariant variant : pmkb.variant()) {
            int idVariant = context.insertInto(PMKBVARIANT,
                    PMKBVARIANT.AMINOACIDCHANGE,
                    PMKBVARIANT.GERMLINE,
                    PMKBVARIANT.PARTNERGENE,
                    PMKBVARIANT.CODONS,
                    PMKBVARIANT.DESCRIPTION,
                    PMKBVARIANT.EXONS,
                    PMKBVARIANT.NOTES,
                    PMKBVARIANT.COSMIC,
                    PMKBVARIANT.EFFECT,
                    PMKBVARIANT.CNVTYPE,
                    PMKBVARIANT.IDVARIANT,
                    PMKBVARIANT.CYTOBAND,
                    PMKBVARIANT.VARIANTTYPE,
                    PMKBVARIANT.DNACHANGE,
                    PMKBVARIANT.COORDINATES,
                    PMKBVARIANT.CHROMOSOMEBASEDCNV,
                    PMKBVARIANT.TRANSCRIPT,
                    PMKBVARIANT.DESCRIPTIONTYPE,
                    PMKBVARIANT.CHROMOSOME,
                    PMKBVARIANT.NAME,
                    PMKBVARIANT.VICCENTRYID)
                    .values(variant.aminoAcidChange(),
                            variant.germline(),
                            variant.partnerGene(),
                            variant.codons(),
                            variant.description(),
                            variant.exons(),
                            variant.notes(),
                            variant.cosmic(),
                            variant.effect(),
                            variant.cnvType(),
                            variant.id(),
                            variant.cytoband(),
                            variant.variantType(),
                            variant.dnaChange(),
                            variant.coordinates(),
                            variant.chromosomeBasedCnv(),
                            variant.transcript(),
                            variant.descriptionType(),
                            variant.chromosome(),
                            variant.name(),
                            viccEntryId)
                    .returning(PMKBVARIANT.ID)
                    .fetchOne()
                    .getValue(PMKBVARIANT.ID);

            for (PmkbGene gene : variant.gene()) {
                context.insertInto(PMKBGENE,
                        PMKBGENE.DESCRIPTION,
                        PMKBGENE.CREATEDAT,
                        PMKBGENE.UPDATEDAT,
                        PMKBGENE.ACTIVEIND,
                        PMKBGENE.EXTERNALID,
                        PMKBGENE.IDGENE,
                        PMKBGENE.NAME,
                        PMKBGENE.PMKBVARIANTID)
                        .values(gene.description(),
                                gene.createdAt(),
                                gene.updatedAt(),
                                gene.activeInd(),
                                gene.externalId(),
                                gene.id(),
                                gene.name(),
                                idVariant)
                        .execute();
            }
        }
    }

    static void deleteAll(@NotNull DSLContext context) {
        context.deleteFrom(PMKB).execute();
        context.deleteFrom(PMKBTISSUE).execute();
        context.deleteFrom(PMKBTUMOR).execute();
        context.deleteFrom(PMKBVARIANT).execute();
        context.deleteFrom(PMKBGENE).execute();
    }
}
