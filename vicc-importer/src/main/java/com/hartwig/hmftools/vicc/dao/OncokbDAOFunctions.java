package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.ONCOKB;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBBIOLOGICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBCLINICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBCONSEQUENCEBIOLOGICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBCONSEQUENCECLINICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBDRUGABSTRACTCLINICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBGENEALIASBIOLOGICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBGENEALIASCLINICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBGENEBIOLOGICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBGENECLINICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBVARIANTBIOLOGICAL;
import static com.hartwig.hmftools.vicc.database.Tables.ONCOKBVARIANTCLINICAL;

import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKb;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKbBiological;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKbClinical;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKbConsequence;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKbDrugAbstract;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKbGene;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKbVariant;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

final class OncokbDAOFunctions {

    private OncokbDAOFunctions() {
    }

    static void write(@NotNull DSLContext context, int viccEntryId, @NotNull OncoKb oncoKb) {
        int id = context.insertInto(ONCOKB, ONCOKB.VICCENTRYID).values(viccEntryId).returning(ONCOKB.ID).fetchOne().getValue(ONCOKB.ID);

        OncoKbBiological oncokbBiological = oncoKb.oncoKbBiological();
        if (oncokbBiological != null) {
            int idBiological = context.insertInto(ONCOKBBIOLOGICAL,
                    ONCOKBBIOLOGICAL.GENE,
                    ONCOKBBIOLOGICAL.ENTREZGENEID,
                    ONCOKBBIOLOGICAL.ISOFORM,
                    ONCOKBBIOLOGICAL.REFSEQ,
                    ONCOKBBIOLOGICAL.ONCOGENIC,
                    ONCOKBBIOLOGICAL.MUTATIONEFFECT,
                    ONCOKBBIOLOGICAL.MUTATIONEFFECTPMIDS,
                    ONCOKBBIOLOGICAL.MUTATIONEFFECTABSTRACTS,
                    ONCOKBBIOLOGICAL.ONCOKBID)
                    .values(oncokbBiological.gene(),
                            oncokbBiological.entrezGeneId(),
                            oncokbBiological.isoform(),
                            oncokbBiological.refSeq(),
                            oncokbBiological.oncogenic(),
                            oncokbBiological.mutationEffect(),
                            oncokbBiological.mutationEffectPmids(),
                            oncokbBiological.mutationEffectAbstracts(),
                            id)
                    .returning(ONCOKBBIOLOGICAL.ID)
                    .fetchOne()
                    .getValue(ONCOKBBIOLOGICAL.ID);

            OncoKbVariant oncokbVariant = oncokbBiological.oncokbVariant();
            int idVariant = context.insertInto(ONCOKBVARIANTBIOLOGICAL,
                    ONCOKBVARIANTBIOLOGICAL.NAME,
                    ONCOKBVARIANTBIOLOGICAL.ALTERATION,
                    ONCOKBVARIANTBIOLOGICAL.PROTEINSTART,
                    ONCOKBVARIANTBIOLOGICAL.PROTEINEND,
                    ONCOKBVARIANTBIOLOGICAL.REFRESIDUES,
                    ONCOKBVARIANTBIOLOGICAL.VARIANTRESIDUES,
                    ONCOKBVARIANTBIOLOGICAL.ONCOKBBIOLOGICALID)
                    .values(oncokbVariant.name(),
                            oncokbVariant.alteration(),
                            oncokbVariant.proteinStart(),
                            oncokbVariant.proteinEnd(),
                            oncokbVariant.refResidues(),
                            oncokbVariant.variantResidues(),
                            idBiological)
                    .returning(ONCOKBVARIANTBIOLOGICAL.ID)
                    .fetchOne()
                    .getValue(ONCOKBVARIANTBIOLOGICAL.ID);

            OncoKbConsequence oncoKbConsequence = oncokbBiological.oncokbVariant().consequence();
            context.insertInto(ONCOKBCONSEQUENCEBIOLOGICAL,
                    ONCOKBCONSEQUENCEBIOLOGICAL.TERM,
                    ONCOKBCONSEQUENCEBIOLOGICAL.DESCRIPTION,
                    ONCOKBCONSEQUENCEBIOLOGICAL.ISGENERALLYTRUNCATING,
                    ONCOKBCONSEQUENCEBIOLOGICAL.ONCOKBVARIANTBIOLOGICALID)
                    .values(oncoKbConsequence.term(), oncoKbConsequence.description(), oncoKbConsequence.isGenerallyTruncating(), idVariant)
                    .execute();

            OncoKbGene oncoKbGene = oncokbBiological.oncokbVariant().gene();
            int idGene = context.insertInto(ONCOKBGENEBIOLOGICAL,
                    ONCOKBGENEBIOLOGICAL.HUGOSYMBOL,
                    ONCOKBGENEBIOLOGICAL.NAME,
                    ONCOKBGENEBIOLOGICAL.ENTREZGENEID,
                    ONCOKBGENEBIOLOGICAL.CURATEDISOFORM,
                    ONCOKBGENEBIOLOGICAL.CURATEDREFSEQ,
                    ONCOKBGENEBIOLOGICAL.ONCOGENE,
                    ONCOKBGENEBIOLOGICAL.TSG,
                    ONCOKBGENEBIOLOGICAL.ONCOKBVARIANTBIOLOGICALID)
                    .values(oncoKbGene.hugoSymbol(),
                            oncoKbGene.name(),
                            oncoKbGene.entrezGeneId(),
                            oncoKbGene.curatedIsoform(),
                            oncoKbGene.curatedRefSeq(),
                            oncoKbGene.oncogene(),
                            oncoKbGene.tsg(),
                            idVariant)
                    .returning(ONCOKBGENEBIOLOGICAL.ID)
                    .fetchOne()
                    .getValue(ONCOKBGENEBIOLOGICAL.ID);

            for (String geneAlias : oncoKbGene.geneAliases()) {
                context.insertInto(ONCOKBGENEALIASBIOLOGICAL,
                        ONCOKBGENEALIASBIOLOGICAL.GENEALIAS,
                        ONCOKBGENEALIASBIOLOGICAL.ONCOKBGENEBIOLOGICALID).values(geneAlias, idGene);
            }
        }

        OncoKbClinical oncokbClinical = oncoKb.oncoKbClinical();
        if (oncokbClinical != null) {
            int idClinical = context.insertInto(ONCOKBCLINICAL,
                    ONCOKBCLINICAL.GENE,
                    ONCOKBCLINICAL.ENTREZGENEID,
                    ONCOKBCLINICAL.ISOFORM,
                    ONCOKBCLINICAL.REFSEQ,
                    ONCOKBCLINICAL.CANCERTYPE,
                    ONCOKBCLINICAL.DRUG,
                    ONCOKBCLINICAL.DRUGPMIDS,
                    ONCOKBCLINICAL.LEVEL,
                    ONCOKBCLINICAL.LEVELLABEL,
                    ONCOKBCLINICAL.ONCOKBID)
                    .values(oncokbClinical.gene(),
                            oncokbClinical.entrezGeneId(),
                            oncokbClinical.isoform(),
                            oncokbClinical.refSeq(),
                            oncokbClinical.cancerType(),
                            oncokbClinical.drug(),
                            oncokbClinical.drugPmids(),
                            oncokbClinical.level(),
                            oncokbClinical.levelLabel(),
                            id)
                    .returning(ONCOKBCLINICAL.ID)
                    .fetchOne()
                    .getValue(ONCOKBCLINICAL.ID);

            for (OncoKbDrugAbstract drugAbstract : oncokbClinical.drugAbstracts()) {
                context.insertInto(ONCOKBDRUGABSTRACTCLINICAL,
                        ONCOKBDRUGABSTRACTCLINICAL.TEXT,
                        ONCOKBDRUGABSTRACTCLINICAL.LINK,
                        ONCOKBDRUGABSTRACTCLINICAL.ONCOKBCLINICALID)
                        .values(drugAbstract.text(), drugAbstract.link(), idClinical)
                        .execute();
            }

            OncoKbVariant variantClinical = oncokbClinical.variant();
            int idClinicalVariant = context.insertInto(ONCOKBVARIANTCLINICAL,
                    ONCOKBVARIANTCLINICAL.NAME,
                    ONCOKBVARIANTCLINICAL.ALTERATION,
                    ONCOKBVARIANTCLINICAL.PROTEINSTART,
                    ONCOKBVARIANTCLINICAL.PROTEINEND,
                    ONCOKBVARIANTCLINICAL.REFRESIDUES,
                    ONCOKBVARIANTCLINICAL.VARIANTRESIDUES,
                    ONCOKBVARIANTCLINICAL.ONCOKBCLINICALID)
                    .values(variantClinical.name(),
                            variantClinical.alteration(),
                            variantClinical.proteinStart(),
                            variantClinical.proteinEnd(),
                            variantClinical.refResidues(),
                            variantClinical.variantResidues(),
                            idClinical)
                    .returning(ONCOKBVARIANTCLINICAL.ID)
                    .fetchOne()
                    .getValue(ONCOKBVARIANTCLINICAL.ID);

            OncoKbConsequence consequenceClinical = oncokbClinical.variant().consequence();
            context.insertInto(ONCOKBCONSEQUENCECLINICAL,
                    ONCOKBCONSEQUENCECLINICAL.TERM,
                    ONCOKBCONSEQUENCECLINICAL.DESCRIPTION,
                    ONCOKBCONSEQUENCECLINICAL.ISGENERALLYTRUNCATING,
                    ONCOKBCONSEQUENCECLINICAL.ONCOKBVARIANTCLINICALID)
                    .values(consequenceClinical.term(),
                            consequenceClinical.description(),
                            consequenceClinical.isGenerallyTruncating(),
                            idClinicalVariant)
                    .execute();

            OncoKbGene geneClinical = oncokbClinical.variant().gene();
            int idGeneClinical = context.insertInto(ONCOKBGENECLINICAL,
                    ONCOKBGENECLINICAL.HUGOSYMBOL,
                    ONCOKBGENECLINICAL.NAME,
                    ONCOKBGENECLINICAL.ENTREZGENEID,
                    ONCOKBGENECLINICAL.CURATEDISOFORM,
                    ONCOKBGENECLINICAL.CURATEDREFSEQ,
                    ONCOKBGENECLINICAL.ONCOGENE,
                    ONCOKBGENECLINICAL.TSG,
                    ONCOKBGENECLINICAL.ONCOKBVARIANTCLINICALID)
                    .values(geneClinical.hugoSymbol(),
                            geneClinical.name(),
                            geneClinical.entrezGeneId(),
                            geneClinical.curatedIsoform(),
                            geneClinical.curatedRefSeq(),
                            geneClinical.oncogene(),
                            geneClinical.tsg(),
                            idClinicalVariant)
                    .returning(ONCOKBGENECLINICAL.ID)
                    .fetchOne()
                    .getValue(ONCOKBGENECLINICAL.ID);

            for (String geneAlias : geneClinical.geneAliases()) {
                context.insertInto(ONCOKBGENEALIASCLINICAL,
                        ONCOKBGENEALIASCLINICAL.GENEALIAS,
                        ONCOKBGENEALIASCLINICAL.ONCOKBGENECLINICALID).values(geneAlias, idGeneClinical).execute();
            }
        }
    }

    static void deleteAll(@NotNull DSLContext context) {
        // First delete biological nodes
        context.deleteFrom(ONCOKBGENEALIASBIOLOGICAL).execute();
        context.deleteFrom(ONCOKBCONSEQUENCEBIOLOGICAL).execute();
        context.deleteFrom(ONCOKBGENEBIOLOGICAL).execute();
        context.deleteFrom(ONCOKBVARIANTBIOLOGICAL).execute();

        // Then delete clinical nodes
        context.deleteFrom(ONCOKBDRUGABSTRACTCLINICAL).execute();
        context.deleteFrom(ONCOKBGENEALIASCLINICAL).execute();
        context.deleteFrom(ONCOKBCONSEQUENCECLINICAL).execute();
        context.deleteFrom(ONCOKBGENECLINICAL).execute();
        context.deleteFrom(ONCOKBVARIANTCLINICAL).execute();

        // Then delete top-nodes + final entry node
        context.deleteFrom(ONCOKBBIOLOGICAL).execute();
        context.deleteFrom(ONCOKBCLINICAL).execute();
        context.deleteFrom(ONCOKB).execute();
    }
}
