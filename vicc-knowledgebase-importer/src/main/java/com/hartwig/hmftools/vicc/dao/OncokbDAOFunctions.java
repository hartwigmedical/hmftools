package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.ONCOKB;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBBIOLOGICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBCLINICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBCONSEQUENCESBIOLOGICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBCONSEQUENCESCLINICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBDRUGABSTRACTSCLINICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBGENEALIASESBIOLOGICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBGENEALIASESCLINICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBGENEBIOLOGICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBGENECLINICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBVARIANTBIOLOGICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBVARIANTCLINICAL;

import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKbBiological;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKbClinical;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKbConsequence;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKbDrugAbstracts;
import com.hartwig.hmftools.vicc.datamodel.oncokb.Oncokb;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncokbGene;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncokbVariant;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

final class OncokbDAOFunctions {

    private OncokbDAOFunctions() {
    }

    static void write(@NotNull DSLContext context, int viccEntryId, @NotNull Oncokb oncokb) {
        int id = context.insertInto(ONCOKB, ONCOKB.VICCENTRYID).values(viccEntryId).returning(ONCOKB.ID).fetchOne().getValue(ONCOKB.ID);

        OncoKbBiological oncokbBiological = oncokb.oncoKbBiological();
        if (oncokbBiological != null) {
            int idBiological = context.insertInto(ONCOKBBIOLOGICAL,
                    ONCOKBBIOLOGICAL.MUTATIONEFFECTPMIDS,
                    ONCOKBBIOLOGICAL.ISOFORM,
                    ONCOKBBIOLOGICAL.ENTREZGENEID,
                    ONCOKBBIOLOGICAL.ONCOGENIC,
                    ONCOKBBIOLOGICAL.MUTATIONEFFECT,
                    ONCOKBBIOLOGICAL.REFSEQ,
                    ONCOKBBIOLOGICAL.GENE,
                    ONCOKBBIOLOGICAL.MUTATIONEFFECTABSTRACTS,
                    ONCOKBBIOLOGICAL.VICCENTRYID)
                    .values(oncokbBiological.mutationEffectPmids(),
                            oncokbBiological.isoform(),
                            oncokbBiological.entrezGeneID(),
                            oncokbBiological.oncogenic(),
                            oncokbBiological.mutationEffect(),
                            oncokbBiological.refSeq(),
                            oncokbBiological.gene(),
                            oncokbBiological.mutationEffectAbstracts(),
                            id)
                    .returning(ONCOKBBIOLOGICAL.ID)
                    .fetchOne()
                    .getValue(ONCOKBBIOLOGICAL.ID);

            OncokbVariant oncokbVariant = oncokbBiological.oncokbVariant();

            int idVariant = context.insertInto(ONCOKBVARIANTBIOLOGICAL,
                    ONCOKBVARIANTBIOLOGICAL.VARIANTRESIDUES,
                    ONCOKBVARIANTBIOLOGICAL.PROTEINSTART,
                    ONCOKBVARIANTBIOLOGICAL.NAME,
                    ONCOKBVARIANTBIOLOGICAL.PROTEINEND,
                    ONCOKBVARIANTBIOLOGICAL.REFRESIDUES,
                    ONCOKBVARIANTBIOLOGICAL.ALTERATION,
                    ONCOKBVARIANTBIOLOGICAL.ONCOKBBIOLOGICALID)
                    .values(oncokbVariant.variantResidues(),
                            oncokbVariant.proteinStart(),
                            oncokbVariant.name(),
                            oncokbVariant.proteinEnd(),
                            oncokbVariant.refResidues(),
                            oncokbVariant.alteration(),
                            idBiological)
                    .returning(ONCOKBVARIANTBIOLOGICAL.ID)
                    .fetchOne()
                    .getValue(ONCOKBVARIANTBIOLOGICAL.ID);

            OncoKbConsequence oncoKbConsequence = oncokbBiological.oncokbVariant().oncoKbConsequence();

            context.insertInto(ONCOKBCONSEQUENCESBIOLOGICAL,
                    ONCOKBCONSEQUENCESBIOLOGICAL.TERM,
                    ONCOKBCONSEQUENCESBIOLOGICAL.DESCRIPTION,
                    ONCOKBCONSEQUENCESBIOLOGICAL.ISGENERALLYTRUNCATING,
                    ONCOKBCONSEQUENCESBIOLOGICAL.ONCOKBVARIANTBIOLOGICALID)
                    .values(oncoKbConsequence.term(), oncoKbConsequence.description(), oncoKbConsequence.isGenerallyTruncating(), idVariant)
                    .execute();

            OncokbGene oncoKbGene = oncokbBiological.oncokbVariant().oncokbGene();

            int idGene = context.insertInto(ONCOKBGENEBIOLOGICAL,
                    ONCOKBGENEBIOLOGICAL.ONCOGENE,
                    ONCOKBGENEBIOLOGICAL.NAME,
                    ONCOKBGENEBIOLOGICAL.HUGOSYMBOL,
                    ONCOKBGENEBIOLOGICAL.CURATEDREFSEQ,
                    ONCOKBGENEBIOLOGICAL.ENTREZGENEID,
                    ONCOKBGENEBIOLOGICAL.TSG,
                    ONCOKBGENEBIOLOGICAL.CURATEDISOFORM,
                    ONCOKBGENEBIOLOGICAL.ONCOKBBIOLOGICALID)
                    .values(oncoKbGene.oncogene(),
                            oncoKbGene.name(),
                            oncoKbGene.hugoSymbol(),
                            oncoKbGene.curatedRefSeq(),
                            oncoKbGene.entrezGeneId(),
                            oncoKbGene.tsg(),
                            oncoKbGene.curatedIsoform(),
                            idVariant)
                    .returning(ONCOKBGENEBIOLOGICAL.ID)
                    .fetchOne()
                    .getValue(ONCOKBGENEBIOLOGICAL.ID);

            for (String geneAliases : oncoKbGene.geneAliases()) {
                context.insertInto(ONCOKBGENEALIASESBIOLOGICAL,
                        ONCOKBGENEALIASESBIOLOGICAL.GENEALIASES,
                        ONCOKBGENEALIASESBIOLOGICAL.ONCOKBGENEBIOLOGICALID).values(geneAliases, idGene);
            }
        }

        OncoKbClinical oncokbClinical = oncokb.oncoKbClinical();
        if (oncokbClinical != null) {
            int idClinical = context.insertInto(ONCOKBCLINICAL,
                    ONCOKBCLINICAL.REFSEQ,
                    ONCOKBCLINICAL.LEVEL,
                    ONCOKBCLINICAL.ENTREZGENEID,
                    ONCOKBCLINICAL.DRUGPMIDS,
                    ONCOKBCLINICAL.CANCERTYPE,
                    ONCOKBCLINICAL.DRUG,
                    ONCOKBCLINICAL.GENE,
                    ONCOKBCLINICAL.LEVELLABEL,
                    ONCOKBCLINICAL.VICCENTRYID)
                    .values(oncokbClinical.refSeq(),
                            oncokbClinical.level(),
                            oncokbClinical.entrezGeneID(),
                            oncokbClinical.drugPmids(),
                            oncokbClinical.cancerType(),
                            oncokbClinical.drug(),
                            oncokbClinical.gene(),
                            oncokbClinical.levelLabel(),
                            id)
                    .returning(ONCOKBCLINICAL.ID)
                    .fetchOne()
                    .getValue(ONCOKBCLINICAL.ID);

            for (OncoKbDrugAbstracts drugAbstracts : oncokbClinical.oncoKbDrugAbstracts()) {
                context.insertInto(ONCOKBDRUGABSTRACTSCLINICAL,
                        ONCOKBDRUGABSTRACTSCLINICAL.TEXT,
                        ONCOKBDRUGABSTRACTSCLINICAL.LINK,
                        ONCOKBDRUGABSTRACTSCLINICAL.ONCOKBCLINICALID)
                        .values(drugAbstracts.text(), drugAbstracts.link(), idClinical)
                        .execute();
            }

            OncokbVariant variantClinical = oncokbClinical.oncokbVariant();
            int idClinicalVariant = context.insertInto(ONCOKBVARIANTCLINICAL,
                    ONCOKBVARIANTCLINICAL.VARIANTRESIDUES,
                    ONCOKBVARIANTCLINICAL.PROTEINSTART,
                    ONCOKBVARIANTCLINICAL.NAME,
                    ONCOKBVARIANTCLINICAL.PROTEINEND,
                    ONCOKBVARIANTCLINICAL.REFRESIDUES,
                    ONCOKBVARIANTCLINICAL.ALTERATION,
                    ONCOKBVARIANTCLINICAL.ONCOKBCLINICALID)
                    .values(variantClinical.variantResidues(),
                            variantClinical.proteinStart(),
                            variantClinical.name(),
                            variantClinical.proteinEnd(),
                            variantClinical.refResidues(),
                            variantClinical.alteration(),
                            idClinical)
                    .returning(ONCOKBVARIANTCLINICAL.ID)
                    .fetchOne()
                    .getValue(ONCOKBVARIANTCLINICAL.ID);

            OncoKbConsequence consequenceClinical = oncokbClinical.oncokbVariant().oncoKbConsequence();

            context.insertInto(ONCOKBCONSEQUENCESCLINICAL,
                    ONCOKBCONSEQUENCESCLINICAL.TERM,
                    ONCOKBCONSEQUENCESCLINICAL.DESCRIPTION,
                    ONCOKBCONSEQUENCESCLINICAL.ISGENERALLYTRUNCATING,
                    ONCOKBCONSEQUENCESCLINICAL.ONCOKBVARIANTCLINICALID)
                    .values(consequenceClinical.term(),
                            consequenceClinical.description(),
                            consequenceClinical.isGenerallyTruncating(),
                            idClinicalVariant)
                    .execute();

            OncokbGene geneClinical = oncokbClinical.oncokbVariant().oncokbGene();

            int idGeneClinical = context.insertInto(ONCOKBGENECLINICAL,
                    ONCOKBGENECLINICAL.ONCOGENE,
                    ONCOKBGENECLINICAL.NAME,
                    ONCOKBGENECLINICAL.HUGOSYMBOL,
                    ONCOKBGENECLINICAL.CURATEDREFSEQ,
                    ONCOKBGENECLINICAL.ENTREZGENEID,
                    ONCOKBGENECLINICAL.TSG,
                    ONCOKBGENECLINICAL.CURATEDISOFORM,
                    ONCOKBGENECLINICAL.ONCOKBCLINICALID)
                    .values(geneClinical.oncogene(),
                            geneClinical.name(),
                            geneClinical.hugoSymbol(),
                            geneClinical.curatedRefSeq(),
                            geneClinical.entrezGeneId(),
                            geneClinical.tsg(),
                            geneClinical.curatedIsoform(),
                            idClinicalVariant)
                    .returning(ONCOKBGENECLINICAL.ID)
                    .fetchOne()
                    .getValue(ONCOKBGENECLINICAL.ID);

            for (String geneAliasesClinical : geneClinical.geneAliases()) {
                context.insertInto(ONCOKBGENEALIASESCLINICAL,
                        ONCOKBGENEALIASESCLINICAL.GENEALIASES,
                        ONCOKBGENEALIASESCLINICAL.ONCOKBGENECLINICALID).values(geneAliasesClinical, idGeneClinical).execute();
            }
        }
    }

    static void deleteAll(@NotNull DSLContext context) {
        // First delete biological nodes
        context.deleteFrom(ONCOKBGENEALIASESBIOLOGICAL).execute();
        context.deleteFrom(ONCOKBCONSEQUENCESBIOLOGICAL).execute();
        context.deleteFrom(ONCOKBGENEBIOLOGICAL).execute();
        context.deleteFrom(ONCOKBVARIANTBIOLOGICAL).execute();

        // Then delete clinical nodes
        context.deleteFrom(ONCOKBDRUGABSTRACTSCLINICAL).execute();
        context.deleteFrom(ONCOKBGENEALIASESCLINICAL).execute();
        context.deleteFrom(ONCOKBCONSEQUENCESCLINICAL).execute();
        context.deleteFrom(ONCOKBGENECLINICAL).execute();
        context.deleteFrom(ONCOKBVARIANTCLINICAL).execute();

        // Then delete top-nodes + final entry node
        context.deleteFrom(ONCOKBBIOLOGICAL).execute();
        context.deleteFrom(ONCOKBCLINICAL).execute();
        context.deleteFrom(ONCOKB).execute();
    }
}
