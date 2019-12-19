package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCH;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHAST;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHASTLEFT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHASTRIGHT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHASTRIGHTLEFT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHASTRIGHTLEFTLEFT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHASTRIGHTLEFTRIGHT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHASTRIGHTRIGHT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONALT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONCHR;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONCOSMICID;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONDBSNP;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONEND;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONEXON;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONEXONICFUNC;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONNUCLEOTIDECHANGE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONPARENTS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONPARENTSTRANSCRIPTS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONPATHOLOGY;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONPOPFREQMAX;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONREF;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONSOURCES;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONSTART;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONTRANSCRIPTS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCRITERIAMET;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCRITERIAUNMET;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHEXTERNALID;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINCLUDECONDITION0;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINCLUDECONDITION1;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINCLUDEDRUG1;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINCLUDEGENE0;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINCLUDEMUTATION0;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINCLUDEMUTATION1;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINCLUDESTAGE0;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINSTUTITION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSCDNA;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON1;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON10;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON11;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON12;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON13;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON14;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON15;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON16;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON17;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON18;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON19;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON2;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON20;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON21;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON22;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON23;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON24;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON25;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON26;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON27;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON28;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON29;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON3;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON30;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON31;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON32;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON33;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON34;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON35;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON36;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON37;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON38;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON39;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON4;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON40;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON41;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON5;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON6;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON7;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON8;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONPOSITIESEXON9;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONSINFO;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIES;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIESEXONPOSITIES;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONABAND;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONACHR;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONACOORD;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONAGENE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONAORI;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONAREG;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONATX;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONBBAND;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONBCHR;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONBCOORD;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONBGENE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONBORI;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONBREG;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONBTX;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSFUSIONINS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSGRCH37LOCATION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSGRCH37LOCATIONCONSEQUENCES;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSGRCH37LOCCONSEQUENCESEXONNUMBER;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSGRCH37LOCCONSEQUENCESTXSITES;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSMUTATIONTYPE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSPARENTS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSPARENTSTRANSCRIPT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSPATHOLOGY;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSSOURCE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSSYNONYMS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCESEXONNUMBER;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSWGSDATA;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSWGSDATACLINVARDBID;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSWGSDATACLINVARDIS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSWGSDATACLINVARSIG;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSWGSDATACLINVARSTATUS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSWGSDATAFULLAA;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSWGSDATAGENE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSWGSMAP;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSWGSMAPPROTCOORDS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSWGSMAPSYNONYMS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHPREFELANCE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHSOURCE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTAGS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTHERAPEUTICCONTEXT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTIEREXPLANATION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHVARIANTINFO;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHVARIANTINFOCONSEQUENCES;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHVARIANTINFOFUSIONS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHVARIANTINFOLOCATIONS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHVARIANTINFOLOCATIONSEXONNUMBER;

import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatch;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAgreg;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchBreg;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchClassification;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchCriteriaUnmet;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchExonsInfo;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchFusion;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchFusionData;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchGRCh37Location;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchLocation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchMutation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchParent;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchPrevalence;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchSource;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchTag;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchTherapeuticContext;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchTierExplanation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchTranscriptConsequence;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchTranscriptConsequencesGRCh37;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchVariantInfo;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchWGSALocation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchWGSAMap;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

final class MolecularMatchDAOFunctions {

    private MolecularMatchDAOFunctions() {
    }

    static void write(@NotNull DSLContext context, int viccEntryId, @NotNull MolecularMatch molecularMatch) {
        // TODO Add gene1, resistance1, finding1 and drugclass1
        int id = context.insertInto(MOLECULARMATCH,
                MOLECULARMATCH.SCORE,
                MOLECULARMATCH.AUTOGENERATENARRATIVE,
                MOLECULARMATCH.CLINICALSIGNIFICANCE,
                MOLECULARMATCH.IDMOLECULARMATCH,
                MOLECULARMATCH.UNIQUEKEY,
                MOLECULARMATCH.CIVICVALUE,
                MOLECULARMATCH.HASHKEY,
                MOLECULARMATCH.REGULATORYBODYAPPROVED,
                MOLECULARMATCH.VERSION,
                MOLECULARMATCH.GUIDELINEBODY,
                MOLECULARMATCH.REGULATORYBODY,
                MOLECULARMATCH.CUSTOMER,
                MOLECULARMATCH.DIRECTION,
                MOLECULARMATCH.AMPCAP,
                MOLECULARMATCH.GUIDELINEVERSION,
                MOLECULARMATCH.TIER,
                MOLECULARMATCH.MVLD,
                MOLECULARMATCH.SIXTIER,
                MOLECULARMATCH.NOTHERAPYAVAILABLE,
                MOLECULARMATCH.NARRATIVE,
                MOLECULARMATCH.EXPRESSION,
                MOLECULARMATCH.BIOMARKERCLASS,
                MOLECULARMATCH.VICCENTRYID)
                .values(molecularMatch.score(),
                        molecularMatch.autoGenerateNarrative(),
                        molecularMatch.clinicalSignificance(),
                        molecularMatch.id(),
                        molecularMatch.uniqueKey(),
                        molecularMatch.civicValue(),
                        molecularMatch.hashKey(),
                        molecularMatch.regulatoryBodyApproved(),
                        molecularMatch.version(),
                        molecularMatch.guidelineBody(),
                        molecularMatch.regulatoryBody(),
                        molecularMatch.customer(),
                        molecularMatch.direction(),
                        molecularMatch.ampcap(),
                        molecularMatch.guidelineVersion(),
                        molecularMatch.tier(),
                        molecularMatch.mvld(),
                        molecularMatch.sixtier(),
                        molecularMatch.noTherapyAvailable(),
                        molecularMatch.narrative(),
                        molecularMatch.expression(),
                        molecularMatch.biomarkerClass(),
                        viccEntryId)
                .returning(MOLECULARMATCH.ID)
                .fetchOne()
                .getValue(MOLECULARMATCH.ID);

        if (molecularMatch.ast() != null) {
            int idAst = context.insertInto(MOLECULARMATCHAST,
                    MOLECULARMATCHAST.RAW,
                    MOLECULARMATCHAST.VALUE,
                    MOLECULARMATCHAST.OPERATOR,
                    MOLECULARMATCHAST.TYPE,
                    MOLECULARMATCHAST.MOLECULARMATCHID)
                    .values(molecularMatch.ast().raw(),
                            molecularMatch.ast().value(),
                            molecularMatch.ast().operator(),
                            molecularMatch.ast().type(),
                            id)
                    .returning(MOLECULARMATCHAST.ID)
                    .fetchOne()
                    .getValue(MOLECULARMATCHAST.ID);

            if (molecularMatch.ast().left() != null) {
                context.insertInto(MOLECULARMATCHASTLEFT,
                        MOLECULARMATCHASTLEFT.RAW,
                        MOLECULARMATCHASTLEFT.VALUE,
                        MOLECULARMATCHASTLEFT.OPERATOR,
                        MOLECULARMATCHASTLEFT.TYPE,
                        MOLECULARMATCHASTLEFT.MOLECULARMATCHASTID)
                        .values(molecularMatch.ast().left().raw(),
                                molecularMatch.ast().left().value(),
                                molecularMatch.ast().left().operator(),
                                molecularMatch.ast().left().type(),
                                idAst)
                        .execute();
            }

            if (molecularMatch.ast().right() != null) {
                int idAstRight = context.insertInto(MOLECULARMATCHASTRIGHT,
                        MOLECULARMATCHASTRIGHT.RAW,
                        MOLECULARMATCHASTRIGHT.VALUE,
                        MOLECULARMATCHASTRIGHT.OPERATOR,
                        MOLECULARMATCHASTRIGHT.TYPE,
                        MOLECULARMATCHASTRIGHT.MOLECULARMATCHASTID)
                        .values(molecularMatch.ast().right().raw(),
                                molecularMatch.ast().right().value(),
                                molecularMatch.ast().right().operator(),
                                molecularMatch.ast().right().type(),
                                idAst)
                        .returning(MOLECULARMATCHASTRIGHT.ID)
                        .fetchOne()
                        .getValue(MOLECULARMATCHASTRIGHT.ID);

                if (molecularMatch.ast().right().right() != null) {
                    context.insertInto(MOLECULARMATCHASTRIGHTRIGHT,
                            MOLECULARMATCHASTRIGHTRIGHT.RAW,
                            MOLECULARMATCHASTRIGHTRIGHT.VALUE,
                            MOLECULARMATCHASTRIGHTRIGHT.TYPE,
                            MOLECULARMATCHASTRIGHTRIGHT.MOLECULARMATCHASTRIGHTID)
                            .values(molecularMatch.ast().right().right().raw(),
                                    molecularMatch.ast().right().right().value(),
                                    molecularMatch.ast().right().right().type(),
                                    idAstRight);
                }

                if (molecularMatch.ast().right().left() != null) {
                    int idRightLeft = context.insertInto(MOLECULARMATCHASTRIGHTLEFT,
                            MOLECULARMATCHASTRIGHTLEFT.RAW,
                            MOLECULARMATCHASTRIGHTLEFT.VALUE,
                            MOLECULARMATCHASTRIGHTLEFT.TYPE,
                            MOLECULARMATCHASTRIGHTLEFT.MOLECULARMATCHASTRIGHTID)
                            .values(molecularMatch.ast().right().left().raw(),
                                    molecularMatch.ast().right().left().value(),
                                    molecularMatch.ast().right().left().type(),
                                    idAstRight)
                            .returning(MOLECULARMATCHASTRIGHTLEFT.ID)
                            .fetchOne()
                            .getValue(MOLECULARMATCHASTRIGHTLEFT.ID);

                    if (molecularMatch.ast().right().left().right() != null) {
                        context.insertInto(MOLECULARMATCHASTRIGHTLEFTRIGHT,
                                MOLECULARMATCHASTRIGHTLEFTRIGHT.RAW,
                                MOLECULARMATCHASTRIGHTLEFTRIGHT.VALUE,
                                MOLECULARMATCHASTRIGHTLEFTRIGHT.TYPE,
                                MOLECULARMATCHASTRIGHTLEFTRIGHT.MOLECULARMATCHASTRIGHTLEFTID)
                                .values(molecularMatch.ast().right().left().right().raw(),
                                        molecularMatch.ast().right().left().right().value(),
                                        molecularMatch.ast().right().left().right().type(),
                                        idRightLeft)
                                .execute();
                    }

                    if (molecularMatch.ast().right().left().left() != null) {
                        context.insertInto(MOLECULARMATCHASTRIGHTLEFTLEFT,
                                MOLECULARMATCHASTRIGHTLEFTLEFT.RAW,
                                MOLECULARMATCHASTRIGHTLEFTLEFT.OPERATOR,
                                MOLECULARMATCHASTRIGHTLEFTLEFT.VALUE,
                                MOLECULARMATCHASTRIGHTLEFTLEFT.TYPE,
                                MOLECULARMATCHASTRIGHTLEFTLEFT.MOLECULARMATCHASTRIGHTLEFTID)
                                .values(molecularMatch.ast().right().left().left().raw(),
                                        molecularMatch.ast().right().left().left().operator(),
                                        molecularMatch.ast().right().left().left().value(),
                                        molecularMatch.ast().right().left().left().type(),
                                        idRightLeft)
                                .execute();
                    }

                }

            }

        }

        if (molecularMatch.institutions() != null) {
            for (String institution : molecularMatch.institutions()) {
                context.insertInto(MOLECULARMATCHINSTUTITION,
                        MOLECULARMATCHINSTUTITION.INSTITUTION,
                        MOLECULARMATCHINSTUTITION.MOLECULARMATCHID).values(institution, id).execute();
            }
        }

        if (molecularMatch.includeGene0() != null) {
            for (String includeGene0 : molecularMatch.includeGene0()) {
                context.insertInto(MOLECULARMATCHINCLUDEGENE0,
                        MOLECULARMATCHINCLUDEGENE0.INCLUDEGENE0,
                        MOLECULARMATCHINCLUDEGENE0.MOLECULARMATCHID).values(includeGene0, id).execute();
            }
        }

        if (molecularMatch.externalIds() != null) {
            for (String externalId : molecularMatch.externalIds()) {
                context.insertInto(MOLECULARMATCHEXTERNALID,
                        MOLECULARMATCHEXTERNALID.EXTERNAL_ID,
                        MOLECULARMATCHEXTERNALID.MOLECULARMATCHID).values(externalId, id).execute();
            }
        }

        if (molecularMatch.includeStage0() != null) {
            for (String includeStage0 : molecularMatch.includeStage0()) {
                context.insertInto(MOLECULARMATCHINCLUDESTAGE0,
                        MOLECULARMATCHINCLUDESTAGE0.INCLUDESTAGE0,
                        MOLECULARMATCHINCLUDESTAGE0.MOLECULARMATCHID).values(includeStage0, id).execute();
            }
        }

        if (molecularMatch.includeDrug1() != null) {
            for (String includeDrug1 : molecularMatch.includeDrug1()) {
                context.insertInto(MOLECULARMATCHINCLUDEDRUG1,
                        MOLECULARMATCHINCLUDEDRUG1.INCLUDEDRUG1,
                        MOLECULARMATCHINCLUDEDRUG1.MOLECULARMATCHID).values(includeDrug1, id).execute();
            }
        }

        for (String includeCondition1 : molecularMatch.includeCondition1()) {
            context.insertInto(MOLECULARMATCHINCLUDECONDITION1,
                    MOLECULARMATCHINCLUDECONDITION1.INCLUDECONDITION1,
                    MOLECULARMATCHINCLUDECONDITION1.MOLECULARMATCHID).values(includeCondition1, id).execute();
        }

        if (molecularMatch.includeMutation1() != null) {
            for (String includeMutation1 : molecularMatch.includeMutation1()) {
                context.insertInto(MOLECULARMATCHINCLUDEMUTATION1,
                        MOLECULARMATCHINCLUDEMUTATION1.INCLUDEMUTATION1,
                        MOLECULARMATCHINCLUDEMUTATION1.MOLECULARMATCHID).values(includeMutation1, id).execute();
            }
        }

        for (String includeCondition0 : molecularMatch.includeCondition0()) {
            context.insertInto(MOLECULARMATCHINCLUDECONDITION0,
                    MOLECULARMATCHINCLUDECONDITION0.INCLUDECONDITION0,
                    MOLECULARMATCHINCLUDECONDITION0.MOLECULARMATCHID).values(includeCondition0, id).execute();
        }

        if (molecularMatch.includeMutation0() != null) {
            for (String includeMutation0 : molecularMatch.includeMutation0()) {
                context.insertInto(MOLECULARMATCHINCLUDEMUTATION0,
                        MOLECULARMATCHINCLUDEMUTATION0.INCLUDEMUTATION0,
                        MOLECULARMATCHINCLUDEMUTATION0.MOLECULARMATCHID).values(includeMutation0, id).execute();
            }
        }

        for (String criteriaMet : molecularMatch.criteriaMets()) {
            context.insertInto(MOLECULARMATCHCRITERIAMET, MOLECULARMATCHCRITERIAMET.CRITERIAMET, MOLECULARMATCHCRITERIAMET.MOLECULARMATCHID)
                    .values(criteriaMet, id)
                    .execute();
        }

        for (MolecularMatchSource source : molecularMatch.sources()) {
            context.insertInto(MOLECULARMATCHSOURCE,
                    MOLECULARMATCHSOURCE.NAME,
                    MOLECULARMATCHSOURCE.SUPPRESS,
                    MOLECULARMATCHSOURCE.PUBID,
                    MOLECULARMATCHSOURCE.SUBTYPE,
                    MOLECULARMATCHSOURCE.VALID,
                    MOLECULARMATCHSOURCE.LINK,
                    MOLECULARMATCHSOURCE.YEAR,
                    MOLECULARMATCHSOURCE.TRIALID,
                    MOLECULARMATCHSOURCE.TYPE,
                    MOLECULARMATCHSOURCE.IDSOURCE,
                    MOLECULARMATCHSOURCE.INSTITUTION,
                    MOLECULARMATCHSOURCE.TRIALPHASE,
                    MOLECULARMATCHSOURCE.FUNCTIONALCONSEQUENCE,
                    MOLECULARMATCHSOURCE.TRUSTRATING,
                    MOLECULARMATCHSOURCE.MOLECULARMATCHID)
                    .values(source.name(),
                            source.suppress(),
                            source.pubId(),
                            source.subType(),
                            source.valid(),
                            source.link(),
                            source.year(),
                            source.trialId(),
                            source.type(),
                            source.id(),
                            source.institution(),
                            source.trialPhase(),
                            source.functionalConsequence(),
                            source.trustRating(),
                            id)
                    .execute();
        }

        for (MolecularMatchTierExplanation tierExplanation : molecularMatch.tierExplanations()) {
            context.insertInto(MOLECULARMATCHTIEREXPLANATION,
                    MOLECULARMATCHTIEREXPLANATION.TIER,
                    MOLECULARMATCHTIEREXPLANATION.STEP,
                    MOLECULARMATCHTIEREXPLANATION.MESSAGE,
                    MOLECULARMATCHTIEREXPLANATION.SUCCESS,
                    MOLECULARMATCHTIEREXPLANATION.MOLECULARMATCHID)
                    .values(tierExplanation.tier(), tierExplanation.step(), tierExplanation.message(), tierExplanation.success(), id)
                    .execute();
        }

        for (MolecularMatchTherapeuticContext therapeuticContext : molecularMatch.therapeuticContexts()) {
            context.insertInto(MOLECULARMATCHTHERAPEUTICCONTEXT,
                    MOLECULARMATCHTHERAPEUTICCONTEXT.FACET,
                    MOLECULARMATCHTHERAPEUTICCONTEXT.NAME,
                    MOLECULARMATCHTHERAPEUTICCONTEXT.SUPPRESS,
                    MOLECULARMATCHTHERAPEUTICCONTEXT.VALID,
                    MOLECULARMATCHTHERAPEUTICCONTEXT.MOLECULARMATCHID)
                    .values(therapeuticContext.facet(),
                            therapeuticContext.name(),
                            therapeuticContext.suppress(),
                            therapeuticContext.valid(),
                            id)
                    .execute();
        }

        for (MolecularMatchTag tags : molecularMatch.tags()) {
            context.insertInto(MOLECULARMATCHTAGS,
                    MOLECULARMATCHTAGS.PRIORITY,
                    MOLECULARMATCHTAGS.COMPOSITEKEY,
                    MOLECULARMATCHTAGS.SUPPRESS,
                    MOLECULARMATCHTAGS.FILTERTYPE,
                    MOLECULARMATCHTAGS.TERM,
                    MOLECULARMATCHTAGS.PRIMARYVALUE,
                    MOLECULARMATCHTAGS.FACET,
                    MOLECULARMATCHTAGS.VALID,
                    MOLECULARMATCHTAGS.CUSTOM,
                    MOLECULARMATCHTAGS.ISNEW,
                    MOLECULARMATCHTAGS.GENERATEDBY,
                    MOLECULARMATCHTAGS.MANUALSUPPRESS,
                    MOLECULARMATCHTAGS.GENERATEDBYTERM,
                    MOLECULARMATCHTAGS.TRANSCRIPT,
                    MOLECULARMATCHTAGS.MOLECULARMATCHID)
                    .values(tags.priority(),
                            tags.compositeKey(),
                            tags.suppress(),
                            tags.filterType(),
                            tags.term(),
                            tags.primary(),
                            tags.facet(),
                            tags.valid(),
                            tags.custom(),
                            tags.isNew(),
                            tags.generatedBy(),
                            tags.manualSuppress(),
                            tags.generatedByTerm(),
                            tags.transcript(),
                            id)
                    .execute();
        }

        for (MolecularMatchVariantInfo variantInfo : molecularMatch.variantInfos()) {
            int idVariantInfo = context.insertInto(MOLECULARMATCHVARIANTINFO,
                    MOLECULARMATCHVARIANTINFO.CLASSIFICATION,
                    MOLECULARMATCHVARIANTINFO.NAME,
                    MOLECULARMATCHVARIANTINFO.GENEFUSIONPARTNER,
                    MOLECULARMATCHVARIANTINFO.COSMIC_ID,
                    MOLECULARMATCHVARIANTINFO.GENE,
                    MOLECULARMATCHVARIANTINFO.TRANSCRIPT,
                    MOLECULARMATCHVARIANTINFO.POPFREQMAX,
                    MOLECULARMATCHVARIANTINFO.MOLECULARMATCHID)
                    .values(variantInfo.classification(),
                            variantInfo.name(),
                            variantInfo.geneFusionPartner(),
                            variantInfo.cosmicId(),
                            variantInfo.gene(),
                            variantInfo.transcript(),
                            variantInfo.popFreqMax(),
                            id)
                    .returning(MOLECULARMATCHVARIANTINFO.ID)
                    .fetchOne()
                    .getValue(MOLECULARMATCHVARIANTINFO.ID);

            for (String consequences : variantInfo.consequences()) {
                context.insertInto(MOLECULARMATCHVARIANTINFOCONSEQUENCES,
                        MOLECULARMATCHVARIANTINFOCONSEQUENCES.CONSEQUENCES,
                        MOLECULARMATCHVARIANTINFOCONSEQUENCES.MOLECULARMATCHVARIANTINFOID).values(consequences, idVariantInfo).execute();
            }

            for (MolecularMatchFusion fusions : variantInfo.fusions()) {
                context.insertInto(MOLECULARMATCHVARIANTINFOFUSIONS,
                        MOLECULARMATCHVARIANTINFOFUSIONS.REFERENCEGENOME,
                        MOLECULARMATCHVARIANTINFOFUSIONS.LBPWREP,
                        MOLECULARMATCHVARIANTINFOFUSIONS.RBPWREP,
                        MOLECULARMATCHVARIANTINFOFUSIONS.EXONNUMBER,
                        MOLECULARMATCHVARIANTINFOFUSIONS.CHR,
                        MOLECULARMATCHVARIANTINFOFUSIONS.RBPWLEP,
                        MOLECULARMATCHVARIANTINFOFUSIONS.INTRONNUMBER,
                        MOLECULARMATCHVARIANTINFOFUSIONS.LBPWLEP,
                        MOLECULARMATCHVARIANTINFOFUSIONS.MOLECULARMATCHVARIANTINFOID)
                        .values(fusions.referenceGenome(),
                                fusions.LBPWREP(),
                                fusions.RBPWREP(),
                                fusions.exonNumber(),
                                fusions.chr(),
                                fusions.RBPWLEP(),
                                fusions.intronNumber(),
                                fusions.RBPWLEP(),
                                idVariantInfo)
                        .execute();
            }

            for (MolecularMatchLocation locations : variantInfo.locations()) {
                int Idlocation = context.insertInto(MOLECULARMATCHVARIANTINFOLOCATIONS,
                        MOLECULARMATCHVARIANTINFOLOCATIONS.AMINOACIDCHANGE,
                        MOLECULARMATCHVARIANTINFOLOCATIONS.INTRONNUMBER,
                        MOLECULARMATCHVARIANTINFOLOCATIONS.STOP,
                        MOLECULARMATCHVARIANTINFOLOCATIONS.START,
                        MOLECULARMATCHVARIANTINFOLOCATIONS.CHR,
                        MOLECULARMATCHVARIANTINFOLOCATIONS.STRAND,
                        MOLECULARMATCHVARIANTINFOLOCATIONS.ALT,
                        MOLECULARMATCHVARIANTINFOLOCATIONS.REFERENCEGENOME,
                        MOLECULARMATCHVARIANTINFOLOCATIONS.REF,
                        MOLECULARMATCHVARIANTINFOLOCATIONS.CDNA,
                        MOLECULARMATCHVARIANTINFOLOCATIONS.MOLECULARMATCHVARIANTINFOID)
                        .values(locations.aminoAcidChange(),
                                locations.intronNumber(),
                                locations.stop(),
                                locations.start(),
                                locations.chr(),
                                locations.strand(),
                                locations.alt(),
                                locations.referenceGenome(),
                                locations.ref(),
                                locations.cdna(),
                                idVariantInfo)
                        .returning(MOLECULARMATCHVARIANTINFOLOCATIONS.ID)
                        .fetchOne()
                        .getValue(MOLECULARMATCHVARIANTINFOLOCATIONS.ID);

                if (locations.exonNumbers() != null) {
                    for (String exonNumber : locations.exonNumbers()) {
                        context.insertInto(MOLECULARMATCHVARIANTINFOLOCATIONSEXONNUMBER,
                                MOLECULARMATCHVARIANTINFOLOCATIONSEXONNUMBER.EXONNUMBER,
                                MOLECULARMATCHVARIANTINFOLOCATIONSEXONNUMBER.MOLECULARMATCHVARIANTINFOLOCATIONSID)
                                .values(exonNumber, Idlocation)
                                .execute();

                    }
                }
            }
        }

        for (MolecularMatchClassification classification : molecularMatch.classifications()) {
            int idClassification = context.insertInto(MOLECULARMATCHCLASSIFICATION,
                    MOLECULARMATCHCLASSIFICATION.CLASSIFICATION,
                    MOLECULARMATCHCLASSIFICATION.CLASSIFICATIONOVERRIDE,
                    MOLECULARMATCHCLASSIFICATION.GENESYMBOL,
                    MOLECULARMATCHCLASSIFICATION.DESCRIPTION,
                    MOLECULARMATCHCLASSIFICATION.PRIORITY,
                    MOLECULARMATCHCLASSIFICATION.EXPANDGENESEARCH,
                    MOLECULARMATCHCLASSIFICATION.DRUGSEXPERIMENTALCOUNT,
                    MOLECULARMATCHCLASSIFICATION.DRUGSAPPROVEDOFFLABELCOUNT,
                    MOLECULARMATCHCLASSIFICATION.COPYNUMBERTYPE,
                    MOLECULARMATCHCLASSIFICATION.PUBLICATIONCOUNT,
                    MOLECULARMATCHCLASSIFICATION.TRANSCRIPT,
                    MOLECULARMATCHCLASSIFICATION.NAME,
                    MOLECULARMATCHCLASSIFICATION.ROOTTERM,
                    MOLECULARMATCHCLASSIFICATION.DRUGSAPPROVEDONLABELCOUNT,
                    MOLECULARMATCHCLASSIFICATION.TRIALCOUNT,
                    MOLECULARMATCHCLASSIFICATION.ALIAS,
                    MOLECULARMATCHCLASSIFICATION.MOLECULARMATCHID)
                    .values(classification.classification(),
                            classification.classificationOverride(),
                            classification.geneSymbol(),
                            classification.description(),
                            classification.priority(),
                            classification.expandGeneSearch(),
                            classification.drugsExperimentalCount(),
                            classification.drugsApprovedOffLabelCount(),
                            classification.copyNumberType(),
                            classification.publicationCount(),
                            classification.transcript(),
                            classification.name(),
                            classification.rootTerm(),
                            classification.drugsApprovedOnLabelCount(),
                            classification.trialCount(),
                            classification.alias(),
                            id)
                    .returning(MOLECULARMATCHCLASSIFICATION.ID)
                    .fetchOne()
                    .getValue(MOLECULARMATCHCLASSIFICATION.ID);

            if (classification.end() != null) {
                for (String end : classification.end()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONEND,
                            MOLECULARMATCHCLASSIFICATIONEND.END,
                            MOLECULARMATCHCLASSIFICATIONEND.MOLECULARMATCHCLASSIFICATIONID).values(end, idClassification).execute();
                }
            }

            if (classification.start() != null) {
                for (String start : classification.start()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONSTART,
                            MOLECULARMATCHCLASSIFICATIONSTART.START,
                            MOLECULARMATCHCLASSIFICATIONSTART.MOLECULARMATCHCLASSIFICATIONID).values(start, idClassification).execute();
                }
            }

            if (classification.chr() != null) {
                for (String chr : classification.chr()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONCHR,
                            MOLECULARMATCHCLASSIFICATIONCHR.CHROMOSOME,
                            MOLECULARMATCHCLASSIFICATIONCHR.MOLECULARMATCHCLASSIFICATIONID).values(chr, idClassification).execute();
                }
            }

            if (classification.pathology() != null) {
                for (String pathology : classification.pathology()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONPATHOLOGY,
                            MOLECULARMATCHCLASSIFICATIONPATHOLOGY.PATHOLOGY,
                            MOLECULARMATCHCLASSIFICATIONPATHOLOGY.MOLECULARMATCHCLASSIFICATIONID)
                            .values(pathology, idClassification)
                            .execute();
                }
            }

            if (classification.ref() != null) {
                for (String ref : classification.ref()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONREF,
                            MOLECULARMATCHCLASSIFICATIONREF.REF,
                            MOLECULARMATCHCLASSIFICATIONREF.MOLECULARMATCHCLASSIFICATIONID).values(ref, idClassification).execute();
                }
            }

            if (classification.exon() != null) {
                for (String exon : classification.exon()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONEXON,
                            MOLECULARMATCHCLASSIFICATIONEXON.EXON,
                            MOLECULARMATCHCLASSIFICATIONEXON.MOLECULARMATCHCLASSIFICATIONID).values(exon, idClassification).execute();
                }
            }

            if (classification.alt() != null) {
                for (String alt : classification.alt()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONALT,
                            MOLECULARMATCHCLASSIFICATIONALT.ALT,
                            MOLECULARMATCHCLASSIFICATIONALT.MOLECULARMATCHCLASSIFICATIONID).values(alt, idClassification).execute();
                }
            }

            if (classification.sources() != null) {
                for (String source : classification.sources()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONSOURCES,
                            MOLECULARMATCHCLASSIFICATIONSOURCES.SOURCE,
                            MOLECULARMATCHCLASSIFICATIONSOURCES.MOLECULARMATCHCLASSIFICATIONID).values(source, idClassification).execute();
                }
            }

            if (classification.transcripts() != null) {
                for (String transcript : classification.transcripts()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONTRANSCRIPTS,
                            MOLECULARMATCHCLASSIFICATIONTRANSCRIPTS.TRANSCRIPT,
                            MOLECULARMATCHCLASSIFICATIONTRANSCRIPTS.MOLECULARMATCHCLASSIFICATIONID)
                            .values(transcript, idClassification)
                            .execute();
                }
            }

            if (classification.cosmicId() != null) {
                for (String cosmicId : classification.cosmicId()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONCOSMICID,
                            MOLECULARMATCHCLASSIFICATIONCOSMICID.COSMIC_ID,
                            MOLECULARMATCHCLASSIFICATIONCOSMICID.MOLECULARMATCHCLASSIFICATIONID)
                            .values(cosmicId, idClassification)
                            .execute();
                }
            }

            if (classification.dbSNP() != null) {
                for (String dbSNP : classification.dbSNP()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONDBSNP,
                            MOLECULARMATCHCLASSIFICATIONDBSNP.DBSNP,
                            MOLECULARMATCHCLASSIFICATIONDBSNP.MOLECULARMATCHCLASSIFICATIONID).values(dbSNP, idClassification).execute();
                }
            }

            if (classification.popFreqMax() != null) {
                for (String popFreqMax : classification.popFreqMax()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONPOPFREQMAX,
                            MOLECULARMATCHCLASSIFICATIONPOPFREQMAX.POPFREQMAX,
                            MOLECULARMATCHCLASSIFICATIONPOPFREQMAX.MOLECULARMATCHCLASSIFICATIONID)
                            .values(popFreqMax, idClassification)
                            .execute();
                }
            }

            if (classification.exonicFunc() != null) {
                for (String exonicFunc : classification.exonicFunc()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONEXONICFUNC,
                            MOLECULARMATCHCLASSIFICATIONEXONICFUNC.EXONICFUNC,
                            MOLECULARMATCHCLASSIFICATIONEXONICFUNC.MOLECULARMATCHCLASSIFICATIONID)
                            .values(exonicFunc, idClassification)
                            .execute();
                }
            }

            if (classification.nucleotideChange() != null) {
                for (String nucleotideChange : classification.nucleotideChange()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONNUCLEOTIDECHANGE,
                            MOLECULARMATCHCLASSIFICATIONNUCLEOTIDECHANGE.NUCLEOTIDECHANGE,
                            MOLECULARMATCHCLASSIFICATIONNUCLEOTIDECHANGE.MOLECULARMATCHCLASSIFICATIONID)
                            .values(nucleotideChange, idClassification)
                            .execute();
                }
            }

            if (classification.parents() != null) {
                for (MolecularMatchParent parents : classification.parents()) {
                    int idParents = context.insertInto(MOLECULARMATCHCLASSIFICATIONPARENTS,
                            MOLECULARMATCHCLASSIFICATIONPARENTS.TYPE,
                            MOLECULARMATCHCLASSIFICATIONPARENTS.NAME,
                            MOLECULARMATCHCLASSIFICATIONPARENTS.ACTIONABLEPARENT,
                            MOLECULARMATCHCLASSIFICATIONPARENTS.MOLECULARMATCHCLASSIFICATIONID)
                            .values(parents.type(), parents.name(), parents.actionableParent(), idClassification)
                            .returning(MOLECULARMATCHCLASSIFICATIONPARENTS.ID)
                            .fetchOne()
                            .getValue(MOLECULARMATCHCLASSIFICATIONPARENTS.ID);

                    for (String transcripts : parents.transcripts()) {
                        context.insertInto(MOLECULARMATCHCLASSIFICATIONPARENTSTRANSCRIPTS,
                                MOLECULARMATCHCLASSIFICATIONPARENTSTRANSCRIPTS.TRANSCRIPT,
                                MOLECULARMATCHCLASSIFICATIONPARENTSTRANSCRIPTS.MOLECULARMATCHCLASSIFICATIONPARENTSID)
                                .values(transcripts, idParents)
                                .execute();
                    }
                }
            }
        }

        for (MolecularMatchCriteriaUnmet criteriaUnmet : molecularMatch.criteriaUnmets()) {
            context.insertInto(MOLECULARMATCHCRITERIAUNMET,
                    MOLECULARMATCHCRITERIAUNMET.PRIORITY,
                    MOLECULARMATCHCRITERIAUNMET.COMPOSITEKEY,
                    MOLECULARMATCHCRITERIAUNMET.ISNEW,
                    MOLECULARMATCHCRITERIAUNMET.GENERATEDBY,
                    MOLECULARMATCHCRITERIAUNMET.MANUALSUPPRESS,
                    MOLECULARMATCHCRITERIAUNMET.GENERATEDBYTERM,
                    MOLECULARMATCHCRITERIAUNMET.SUPPRESS,
                    MOLECULARMATCHCRITERIAUNMET.FILTERTYPE,
                    MOLECULARMATCHCRITERIAUNMET.TERM,
                    MOLECULARMATCHCRITERIAUNMET.PRIMARYVALUE,
                    MOLECULARMATCHCRITERIAUNMET.FACET,
                    MOLECULARMATCHCRITERIAUNMET.VALID,
                    MOLECULARMATCHCRITERIAUNMET.CUSTOM,
                    MOLECULARMATCHCRITERIAUNMET.TRANSCRIPT,
                    MOLECULARMATCHCRITERIAUNMET.MOLECULARMATCHID)
                    .values(criteriaUnmet.priority(),
                            criteriaUnmet.compositeKey(),
                            criteriaUnmet.isNew(),
                            criteriaUnmet.generatedBy(),
                            criteriaUnmet.manualSuppress(),
                            criteriaUnmet.generatedByTerm(),
                            criteriaUnmet.suppress(),
                            criteriaUnmet.filterType(),
                            criteriaUnmet.term(),
                            criteriaUnmet.primary(),
                            criteriaUnmet.facet(),
                            criteriaUnmet.valid(),
                            criteriaUnmet.custom(),
                            criteriaUnmet.transcript(),
                            id)
                    .execute();
        }

        for (MolecularMatchPrevalence prevelance : molecularMatch.prevalences()) {
            context.insertInto(MOLECULARMATCHPREFELANCE,
                    MOLECULARMATCHPREFELANCE.COUNT,
                    MOLECULARMATCHPREFELANCE.PERCENT,
                    MOLECULARMATCHPREFELANCE.STUDYID,
                    MOLECULARMATCHPREFELANCE.SAMPLES,
                    MOLECULARMATCHPREFELANCE.MOLECULAR,
                    MOLECULARMATCHPREFELANCE.CONDITIONVALUE,
                    MOLECULARMATCHPREFELANCE.MOLECULARMATCHID)
                    .values(prevelance.count(),
                            prevelance.percent(),
                            prevelance.studyId(),
                            prevelance.samples(),
                            prevelance.molecular(),
                            prevelance.condition(),
                            id)
                    .execute();
        }

        for (MolecularMatchMutation mutations : molecularMatch.mutations()) {
            int idMutations = context.insertInto(MOLECULARMATCHMUTATIONS,
                    MOLECULARMATCHMUTATIONS.LONGESTTRANSCRIPT,
                    MOLECULARMATCHMUTATIONS.DESCRIPTION,
                    MOLECULARMATCHMUTATIONS.SRC,
                    MOLECULARMATCHMUTATIONS.UNIPROTTRANSCRIPT,
                    MOLECULARMATCHMUTATIONS.TRANSCRIPTRECOGNIZED,
                    MOLECULARMATCHMUTATIONS.GENESYMBOL,
                    MOLECULARMATCHMUTATIONS.TRANSCRIPT,
                    MOLECULARMATCHMUTATIONS.IDMUTATIONS,
                    MOLECULARMATCHMUTATIONS.NAME,
                    MOLECULARMATCHMUTATIONS.MOLECULARMATCHID)
                    .values(mutations.longestTranscript(),
                            mutations.description(),
                            mutations.src(),
                            mutations.uniprotTranscript(),
                            mutations.transcriptRecognized(),
                            mutations.geneSymbol(),
                            mutations.transcript(),
                            mutations.id(),
                            mutations.name(),
                            id)
                    .returning(MOLECULARMATCHMUTATIONS.ID)
                    .fetchOne()
                    .getValue(MOLECULARMATCHMUTATIONS.ID);

            for (String mutationType : mutations.mutationTypes()) {
                context.insertInto(MOLECULARMATCHMUTATIONSMUTATIONTYPE,
                        MOLECULARMATCHMUTATIONSMUTATIONTYPE.MUTATIONTYPE,
                        MOLECULARMATCHMUTATIONSMUTATIONTYPE.MOLECULARMATCHMUTATIONSID).values(mutationType, idMutations).execute();
            }

            for (String source : mutations.sources()) {
                context.insertInto(MOLECULARMATCHMUTATIONSSOURCE,
                        MOLECULARMATCHMUTATIONSSOURCE.SOURCE,
                        MOLECULARMATCHMUTATIONSSOURCE.MOLECULARMATCHMUTATIONSID).values(source, idMutations).execute();
            }

            for (String synonyms : mutations.synonyms()) {
                context.insertInto(MOLECULARMATCHMUTATIONSSYNONYMS,
                        MOLECULARMATCHMUTATIONSSYNONYMS.SYNONYMS,
                        MOLECULARMATCHMUTATIONSSYNONYMS.MOLECULARMATCHMUTATIONSID).values(synonyms, idMutations).execute();
            }

            for (String pathology : mutations.pathology()) {
                context.insertInto(MOLECULARMATCHMUTATIONSPATHOLOGY,
                        MOLECULARMATCHMUTATIONSPATHOLOGY.PATHOLOGY,
                        MOLECULARMATCHMUTATIONSPATHOLOGY.MOLECULARMATCHMUTATIONSID).values(pathology, idMutations).execute();
            }

            for (String cDNA : mutations.cDNA()) {
                context.insertInto(MOLECULARMATCHMUTATIONSCDNA,
                        MOLECULARMATCHMUTATIONSCDNA.CDNA,
                        MOLECULARMATCHMUTATIONSCDNA.MOLECULARMATCHMUTATIONSID).values(cDNA, idMutations).execute();
            }

            for (MolecularMatchTranscriptConsequence transcriptConsequence : mutations.transcriptConsequences()) {
                int idTranscriptConsequence = context.insertInto(MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.AMINOACIDCHANGE,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.COMPOSITEKEY,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.INTRONNUMBER,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.SUPPRESS,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.STOP,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.CUSTOM,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.START,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.CHR,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.STRAND,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.VALIDATED,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.TRANSCRIPT,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.CDNA,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.REFERENCEGENOME,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.REF,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.ALT,
                        MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.MOLECULARMATCHMUTATIONSID)
                        .values(transcriptConsequence.aminoAcidChange(),
                                transcriptConsequence.compositeKey(),
                                transcriptConsequence.intronNumber(),
                                transcriptConsequence.suppress(),
                                transcriptConsequence.stop(),
                                transcriptConsequence.custom(),
                                transcriptConsequence.start(),
                                transcriptConsequence.chr(),
                                transcriptConsequence.strand(),
                                transcriptConsequence.validated(),
                                transcriptConsequence.transcript(),
                                transcriptConsequence.cdna(),
                                transcriptConsequence.referenceGenome(),
                                transcriptConsequence.ref(),
                                transcriptConsequence.alt(),
                                idMutations)
                        .returning(MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.ID)
                        .fetchOne()
                        .getValue(MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES.ID);

                if (transcriptConsequence.exonNumbers() != null) {
                    for (String exonNumber : transcriptConsequence.exonNumbers()) {
                        context.insertInto(MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCESEXONNUMBER,
                                MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCESEXONNUMBER.EXONNUMBER,
                                MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCESEXONNUMBER.MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCESID)
                                .values(exonNumber, idTranscriptConsequence)
                                .execute();
                    }
                }
            }

            for (MolecularMatchParent parents : mutations.parents()) {
                int idParents = context.insertInto(MOLECULARMATCHMUTATIONSPARENTS,
                        MOLECULARMATCHMUTATIONSPARENTS.TYPE,
                        MOLECULARMATCHMUTATIONSPARENTS.NAME,
                        MOLECULARMATCHMUTATIONSPARENTS.ACTIONABLEPARENT,
                        MOLECULARMATCHMUTATIONSPARENTS.MOLECULARMATCHMUTATIONSID)
                        .values(parents.type(), parents.name(), parents.actionableParent(), idMutations)
                        .returning(MOLECULARMATCHMUTATIONSPARENTS.ID)
                        .fetchOne()
                        .getValue(MOLECULARMATCHMUTATIONSPARENTS.ID);

                for (String transcripts : parents.transcripts()) {
                    context.insertInto(MOLECULARMATCHMUTATIONSPARENTSTRANSCRIPT,
                            MOLECULARMATCHMUTATIONSPARENTSTRANSCRIPT.TRANSCRIPT,
                            MOLECULARMATCHMUTATIONSPARENTSTRANSCRIPT.MOLECULARMATCHMUTATIONSPARENTSID)
                            .values(transcripts, idParents)
                            .execute();
                }
            }

            if (mutations.wgsaLocations() != null) {
                for (MolecularMatchWGSALocation wgsaDataLocation : mutations.wgsaLocations()) {
                    int idWGSData = context.insertInto(MOLECULARMATCHMUTATIONSWGSDATA,
                            MOLECULARMATCHMUTATIONSWGSDATA.EXONICFUNC,
                            MOLECULARMATCHMUTATIONSWGSDATA.DBSNP,
                            MOLECULARMATCHMUTATIONSWGSDATA.EXAC_NFE,
                            MOLECULARMATCHMUTATIONSWGSDATA.EXAC_FIN,
                            MOLECULARMATCHMUTATIONSWGSDATA.G1000_ALL,
                            MOLECULARMATCHMUTATIONSWGSDATA.G1000_SAS,
                            MOLECULARMATCHMUTATIONSWGSDATA.G1000_EAS,
                            MOLECULARMATCHMUTATIONSWGSDATA.G1000_AFR,
                            MOLECULARMATCHMUTATIONSWGSDATA.EXAC_SAS,
                            MOLECULARMATCHMUTATIONSWGSDATA.EXAC_EAS,
                            MOLECULARMATCHMUTATIONSWGSDATA.EXAC_AMR,
                            MOLECULARMATCHMUTATIONSWGSDATA.EXAC_AFR,
                            MOLECULARMATCHMUTATIONSWGSDATA.EXAC_FREQ,
                            MOLECULARMATCHMUTATIONSWGSDATA.END,
                            MOLECULARMATCHMUTATIONSWGSDATA.START,
                            MOLECULARMATCHMUTATIONSWGSDATA.SIPHY_29WAY_LOGODDS,
                            MOLECULARMATCHMUTATIONSWGSDATA.REF,
                            MOLECULARMATCHMUTATIONSWGSDATA.GERP_RS,
                            MOLECULARMATCHMUTATIONSWGSDATA.FATHMM,
                            MOLECULARMATCHMUTATIONSWGSDATA.NUCLEOTIDECHANGE,
                            MOLECULARMATCHMUTATIONSWGSDATA.PHYLOP100WAY_VERTEBRATE,
                            MOLECULARMATCHMUTATIONSWGSDATA.FUNC,
                            MOLECULARMATCHMUTATIONSWGSDATA.GWAS_PUBMED,
                            MOLECULARMATCHMUTATIONSWGSDATA.TRANSCRIPT,
                            MOLECULARMATCHMUTATIONSWGSDATA.ESP6500SI_AA,
                            MOLECULARMATCHMUTATIONSWGSDATA.ESP6500SI_EA,
                            MOLECULARMATCHMUTATIONSWGSDATA.G1000_EUR,
                            MOLECULARMATCHMUTATIONSWGSDATA.G1000_AMR,
                            MOLECULARMATCHMUTATIONSWGSDATA.CHR_START_REF_ALT,
                            MOLECULARMATCHMUTATIONSWGSDATA.AA,
                            MOLECULARMATCHMUTATIONSWGSDATA.POPFREQMAX,
                            MOLECULARMATCHMUTATIONSWGSDATA.FATHMM_PRED,
                            MOLECULARMATCHMUTATIONSWGSDATA.WGRNA,
                            MOLECULARMATCHMUTATIONSWGSDATA.PHYLOP46WAY_PLACENTAL,
                            MOLECULARMATCHMUTATIONSWGSDATA.KEYVALUE,
                            MOLECULARMATCHMUTATIONSWGSDATA.TARGETSCANS,
                            MOLECULARMATCHMUTATIONSWGSDATA.CHR,
                            MOLECULARMATCHMUTATIONSWGSDATA.COSMIC_ID,
                            MOLECULARMATCHMUTATIONSWGSDATA.ALT,
                            MOLECULARMATCHMUTATIONSWGSDATA.GWAS_DIS,
                            MOLECULARMATCHMUTATIONSWGSDATA.GWAS_SNP,
                            MOLECULARMATCHMUTATIONSWGSDATA.MOLECULARMATCHMUTATIONSID)
                            .values(wgsaDataLocation.exonicFunc(),
                                    wgsaDataLocation.dbSNP(),
                                    wgsaDataLocation.exacNFE(),
                                    wgsaDataLocation.exacFIN(),
                                    wgsaDataLocation.g1000ALL(),
                                    wgsaDataLocation.g1000SAS(),
                                    wgsaDataLocation.g1000EAS(),
                                    wgsaDataLocation.g1000AFR(),
                                    wgsaDataLocation.exacSAS(),
                                    wgsaDataLocation.exacEAS(),
                                    wgsaDataLocation.exacAMR(),
                                    wgsaDataLocation.exacAFR(),
                                    wgsaDataLocation.exacFreq(),
                                    wgsaDataLocation.end(),
                                    wgsaDataLocation.start(),
                                    wgsaDataLocation.siPhy29wayLogOdds(),
                                    wgsaDataLocation.ref(),
                                    wgsaDataLocation.gerpRS(),
                                    wgsaDataLocation.fathmm(),
                                    wgsaDataLocation.nucleotideChange(),
                                    wgsaDataLocation.phyloP100wayVertebrate(),
                                    wgsaDataLocation.func(),
                                    wgsaDataLocation.gwasPubmed(),
                                    wgsaDataLocation.transcript(),
                                    wgsaDataLocation.esp6500siAA(),
                                    wgsaDataLocation.esp6500siEA(),
                                    wgsaDataLocation.g1000EUR(),
                                    wgsaDataLocation.g1000AMR(),
                                    wgsaDataLocation.chrStartRefAlt(),
                                    wgsaDataLocation.aa(),
                                    wgsaDataLocation.popFreqMax(),
                                    wgsaDataLocation.fathmmPred(),
                                    wgsaDataLocation.wgRna(),
                                    wgsaDataLocation.phyloP46wayPlacental(),
                                    wgsaDataLocation.key(),
                                    wgsaDataLocation.targetScanS(),
                                    wgsaDataLocation.chr(),
                                    wgsaDataLocation.cosmicId(),
                                    wgsaDataLocation.alt(),
                                    wgsaDataLocation.gwasDIS(),
                                    wgsaDataLocation.gwasSNP(),
                                    idMutations)
                            .returning(MOLECULARMATCHMUTATIONSWGSDATA.ID)
                            .fetchOne()
                            .getValue(MOLECULARMATCHMUTATIONSWGSDATA.ID);

                    if (wgsaDataLocation.clinVarDiseases() != null) {
                        for (String clinvarDIS : wgsaDataLocation.clinVarDiseases()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSWGSDATACLINVARDIS,
                                    MOLECULARMATCHMUTATIONSWGSDATACLINVARDIS.CLINVARDIS,
                                    MOLECULARMATCHMUTATIONSWGSDATACLINVARDIS.MOLECULARMATCHMUTATIONSWGSDATAID)
                                    .values(clinvarDIS, idWGSData)
                                    .execute();
                        }
                    }

                    if (wgsaDataLocation.clinVarSigs() != null) {
                        for (String clinvarSig : wgsaDataLocation.clinVarSigs()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSWGSDATACLINVARSIG,
                                    MOLECULARMATCHMUTATIONSWGSDATACLINVARSIG.CLINVARSIG,
                                    MOLECULARMATCHMUTATIONSWGSDATACLINVARSIG.MOLECULARMATCHMUTATIONSWGSDATAID)
                                    .values(clinvarSig, idWGSData)
                                    .execute();
                        }
                    }

                    if (wgsaDataLocation.clinVarStates() != null) {
                        for (String clinvarStatus : wgsaDataLocation.clinVarStates()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSWGSDATACLINVARSTATUS,
                                    MOLECULARMATCHMUTATIONSWGSDATACLINVARSTATUS.CLINVARSTATUS,
                                    MOLECULARMATCHMUTATIONSWGSDATACLINVARSTATUS.MOLECULARMATCHMUTATIONSWGSDATAID)
                                    .values(clinvarStatus, idWGSData)
                                    .execute();
                        }
                    }

                    if (wgsaDataLocation.clinVarDbIds() != null) {
                        for (String clinvarDBID : wgsaDataLocation.clinVarDbIds()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSWGSDATACLINVARDBID,
                                    MOLECULARMATCHMUTATIONSWGSDATACLINVARDBID.CLINVARDBID,
                                    MOLECULARMATCHMUTATIONSWGSDATACLINVARDBID.MOLECULARMATCHMUTATIONSWGSDATAID)
                                    .values(clinvarDBID, idWGSData)
                                    .execute();
                        }
                    }

                    for (String fullAA : wgsaDataLocation.fullAAs()) {
                        context.insertInto(MOLECULARMATCHMUTATIONSWGSDATAFULLAA,
                                MOLECULARMATCHMUTATIONSWGSDATAFULLAA.MOLECULARMATCHMUTATIONSWGSDATAFULLAA.FULLAA,
                                MOLECULARMATCHMUTATIONSWGSDATAFULLAA.MOLECULARMATCHMUTATIONSWGSDATAID).values(fullAA, idWGSData).execute();
                    }

                    for (String gene : wgsaDataLocation.genes()) {
                        context.insertInto(MOLECULARMATCHMUTATIONSWGSDATAGENE,
                                MOLECULARMATCHMUTATIONSWGSDATAGENE.GENE,
                                MOLECULARMATCHMUTATIONSWGSDATAGENE.MOLECULARMATCHMUTATIONSWGSDATAID).values(gene, idWGSData).execute();
                    }
                }
            }

            if (mutations.wgsaMaps() != null) {
                for (MolecularMatchWGSAMap wgSaMap : mutations.wgsaMaps()) {
                    int idMap = context.insertInto(MOLECULARMATCHMUTATIONSWGSMAP,
                            MOLECULARMATCHMUTATIONSWGSMAP.AA,
                            MOLECULARMATCHMUTATIONSWGSMAP.NAME,
                            MOLECULARMATCHMUTATIONSWGSMAP.GRCH37_CHR_START_REF_ALT,
                            MOLECULARMATCHMUTATIONSWGSMAP.NUCLEOTIDECHANGE,
                            MOLECULARMATCHMUTATIONSWGSMAP.EXON,
                            MOLECULARMATCHMUTATIONSWGSMAP.GENE,
                            MOLECULARMATCHMUTATIONSWGSMAP.TRANSCRIPT,
                            MOLECULARMATCHMUTATIONSWGSMAP.MOLECULARMATCHMUTATIONSID)
                            .values(wgSaMap.aa(),
                                    wgSaMap.name(),
                                    wgSaMap.grch37ChrStartRefAlt(),
                                    wgSaMap.nucleotideChange(),
                                    wgSaMap.exon(),
                                    wgSaMap.gene(),
                                    wgSaMap.transcript(),
                                    idMutations)
                            .returning(MOLECULARMATCHMUTATIONSWGSMAP.ID)
                            .fetchOne()
                            .getValue(MOLECULARMATCHMUTATIONSWGSMAP.ID);

                    for (String synonyms : wgSaMap.synonyms()) {
                        context.insertInto(MOLECULARMATCHMUTATIONSWGSMAPSYNONYMS,
                                MOLECULARMATCHMUTATIONSWGSMAPSYNONYMS.SYNONYMS,
                                MOLECULARMATCHMUTATIONSWGSMAPSYNONYMS.MOLECULARMATCHMUTATIONSWGSMAPID).values(synonyms, idMap).execute();
                    }

                    for (String protCoords : wgSaMap.protCoords()) {
                        context.insertInto(MOLECULARMATCHMUTATIONSWGSMAPPROTCOORDS,
                                MOLECULARMATCHMUTATIONSWGSMAPPROTCOORDS.PROTCOORDS,
                                MOLECULARMATCHMUTATIONSWGSMAPPROTCOORDS.MOLECULARMATCHMUTATIONSWGSMAPID)
                                .values(protCoords, idMap)
                                .execute();
                    }

                }
            }

            for (MolecularMatchGRCh37Location gRch37Location : mutations.grch37Locations()) {
                int idLocation = context.insertInto(MOLECULARMATCHMUTATIONSGRCH37LOCATION,
                        MOLECULARMATCHMUTATIONSGRCH37LOCATION.COMPOSITEKEY,
                        MOLECULARMATCHMUTATIONSGRCH37LOCATION.REF,
                        MOLECULARMATCHMUTATIONSGRCH37LOCATION.STOP,
                        MOLECULARMATCHMUTATIONSGRCH37LOCATION.START,
                        MOLECULARMATCHMUTATIONSGRCH37LOCATION.CHR,
                        MOLECULARMATCHMUTATIONSGRCH37LOCATION.ALT,
                        MOLECULARMATCHMUTATIONSGRCH37LOCATION.VALIDATED,
                        MOLECULARMATCHMUTATIONSGRCH37LOCATION.STRAND,
                        MOLECULARMATCHMUTATIONSGRCH37LOCATION.MOLECULARMATCHMUTATIONSID)
                        .values(gRch37Location.compositeKey(),
                                gRch37Location.ref(),
                                gRch37Location.stop(),
                                gRch37Location.start(),
                                gRch37Location.chr(),
                                gRch37Location.alt(),
                                gRch37Location.validated(),
                                gRch37Location.strand(),
                                idMutations)
                        .returning(MOLECULARMATCHMUTATIONSGRCH37LOCATION.ID)
                        .fetchOne()
                        .getValue(MOLECULARMATCHMUTATIONSGRCH37LOCATION.ID);

                for (MolecularMatchTranscriptConsequencesGRCh37 transcriptConsequencesGRCH37 : gRch37Location.transcriptConsequences()) {
                    int idConsequences = context.insertInto(MOLECULARMATCHMUTATIONSGRCH37LOCATIONCONSEQUENCES,
                            MOLECULARMATCHMUTATIONSGRCH37LOCATIONCONSEQUENCES.AMINOACIDCHANGE,
                            MOLECULARMATCHMUTATIONSGRCH37LOCATIONCONSEQUENCES.INTRONNUMBER,
                            MOLECULARMATCHMUTATIONSGRCH37LOCATIONCONSEQUENCES.TRANSCRIPT,
                            MOLECULARMATCHMUTATIONSGRCH37LOCATIONCONSEQUENCES.CDNA,
                            MOLECULARMATCHMUTATIONSGRCH37LOCATIONCONSEQUENCES.MOLECULARMATCHMUTATIONSGRCH37LOCATIONID)
                            .values(transcriptConsequencesGRCH37.aminoAcidChange(),
                                    transcriptConsequencesGRCH37.intronNumber(),
                                    transcriptConsequencesGRCH37.transcript(),
                                    transcriptConsequencesGRCH37.cdna(),
                                    idLocation)
                            .returning(MOLECULARMATCHMUTATIONSGRCH37LOCATIONCONSEQUENCES.ID)
                            .fetchOne()
                            .getValue(MOLECULARMATCHMUTATIONSGRCH37LOCATIONCONSEQUENCES.ID);

                    for (String txSites : transcriptConsequencesGRCH37.txSites()) {
                        context.insertInto(MOLECULARMATCHMUTATIONSGRCH37LOCCONSEQUENCESTXSITES,
                                MOLECULARMATCHMUTATIONSGRCH37LOCCONSEQUENCESTXSITES.TXSITES,
                                MOLECULARMATCHMUTATIONSGRCH37LOCCONSEQUENCESTXSITES.MOLECULARMATCHMUTATIONSGRCH37LOCATIONCONSEQUENCESID)
                                .values(txSites, idConsequences)
                                .execute();
                    }

                    if (transcriptConsequencesGRCH37.exonNumbers() != null) {
                        for (String exonNumber : transcriptConsequencesGRCH37.exonNumbers()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSGRCH37LOCCONSEQUENCESEXONNUMBER,
                                    MOLECULARMATCHMUTATIONSGRCH37LOCCONSEQUENCESEXONNUMBER.EXONNUMBER,
                                    MOLECULARMATCHMUTATIONSGRCH37LOCCONSEQUENCESEXONNUMBER.MOLECULARMATCHMUTATIONSGRCH37LOCATIONCONSEQUENCESID)
                                    .values(exonNumber, idConsequences)
                                    .execute();
                        }
                    }
                }
            }

            if (mutations.fusionData() != null) {
                for (MolecularMatchFusionData fusions : mutations.fusionData()) {
                    int idFusions = context.insertInto(MOLECULARMATCHMUTATIONSFUSION,
                            MOLECULARMATCHMUTATIONSFUSION.SYNONYM,
                            MOLECULARMATCHMUTATIONSFUSION.SOURCE,
                            MOLECULARMATCHMUTATIONSFUSION.PAPER,
                            MOLECULARMATCHMUTATIONSFUSION.MOLECULARMATCHMUTATIONSID)
                            .values(fusions.synonym(), fusions.source(), fusions.Paper(), idMutations)
                            .returning(MOLECULARMATCHMUTATIONSFUSION.ID)
                            .fetchOne()
                            .getValue(MOLECULARMATCHMUTATIONSFUSION.ID);

                    if (fusions.Bchr() != null) {
                        for (String Bchr : fusions.Bchr()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONBCHR,
                                    MOLECULARMATCHMUTATIONSFUSIONBCHR.BRCHR,
                                    MOLECULARMATCHMUTATIONSFUSIONBCHR.MOLECULARMATCHMUTATIONSFUSIONID).values(Bchr, idFusions).execute();
                        }
                    }

                    if (fusions.Agene() != null) {
                        for (String agene : fusions.Agene()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONAGENE,
                                    MOLECULARMATCHMUTATIONSFUSIONAGENE.AGENE,
                                    MOLECULARMATCHMUTATIONSFUSIONAGENE.MOLECULARMATCHMUTATIONSFUSIONID).values(agene, idFusions).execute();
                        }
                    }

                    if (fusions.Btx() != null) {
                        for (String btx : fusions.Btx()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONBTX,
                                    MOLECULARMATCHMUTATIONSFUSIONBTX.BTX,
                                    MOLECULARMATCHMUTATIONSFUSIONBTX.MOLECULARMATCHMUTATIONSFUSIONID).values(btx, idFusions).execute();
                        }
                    }

                    if (fusions.Achr() != null) {
                        for (String achr : fusions.Achr()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONATX,
                                    MOLECULARMATCHMUTATIONSFUSIONATX.ATX,
                                    MOLECULARMATCHMUTATIONSFUSIONATX.MOLECULARMATCHMUTATIONSFUSIONID).values(achr, idFusions).execute();
                        }
                    }

                    if (fusions.ins() != null) {
                        for (String ins : fusions.ins()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONINS,
                                    MOLECULARMATCHMUTATIONSFUSIONINS.INS,
                                    MOLECULARMATCHMUTATIONSFUSIONINS.MOLECULARMATCHMUTATIONSFUSIONID).values(ins, idFusions).execute();
                        }
                    }

                    if (fusions.Bgene() != null) {
                        for (String bgene : fusions.Bgene()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONBGENE,
                                    MOLECULARMATCHMUTATIONSFUSIONBGENE.BGENE,
                                    MOLECULARMATCHMUTATIONSFUSIONBGENE.MOLECULARMATCHMUTATIONSFUSIONID).values(bgene, idFusions).execute();
                        }
                    }

                    if (fusions.Acoord() != null) {
                        for (String Acoord : fusions.Acoord()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONACOORD,
                                    MOLECULARMATCHMUTATIONSFUSIONACOORD.ACOORD,
                                    MOLECULARMATCHMUTATIONSFUSIONACOORD.MOLECULARMATCHMUTATIONSFUSIONID)
                                    .values(Acoord, idFusions)
                                    .execute();
                        }
                    }

                    if (fusions.Bori() != null) {
                        for (String bori : fusions.Bori()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONBORI,
                                    MOLECULARMATCHMUTATIONSFUSIONBORI.BORI,
                                    MOLECULARMATCHMUTATIONSFUSIONBORI.MOLECULARMATCHMUTATIONSFUSIONID).values(bori, idFusions).execute();
                        }
                    }

                    if (fusions.Aband() != null) {
                        for (String aband : fusions.Aband()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONABAND,
                                    MOLECULARMATCHMUTATIONSFUSIONABAND.ABAND,
                                    MOLECULARMATCHMUTATIONSFUSIONABAND.MOLECULARMATCHMUTATIONSFUSIONID).values(aband, idFusions).execute();
                        }
                    }

                    if (fusions.Bband() != null) {
                        for (String bband : fusions.Bband()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONBBAND,
                                    MOLECULARMATCHMUTATIONSFUSIONBBAND.BBAND,
                                    MOLECULARMATCHMUTATIONSFUSIONBBAND.MOLECULARMATCHMUTATIONSFUSIONID).values(bband, idFusions).execute();
                        }
                    }

                    if (fusions.Aori() != null) {
                        for (String aori : fusions.Aori()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONAORI,
                                    MOLECULARMATCHMUTATIONSFUSIONAORI.AORI,
                                    MOLECULARMATCHMUTATIONSFUSIONAORI.MOLECULARMATCHMUTATIONSFUSIONID).values(aori, idFusions).execute();
                        }
                    }

                    if (fusions.Atx() != null) {
                        for (String atx : fusions.Atx()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONATX,
                                    MOLECULARMATCHMUTATIONSFUSIONATX.ATX,
                                    MOLECULARMATCHMUTATIONSFUSIONATX.MOLECULARMATCHMUTATIONSFUSIONID).values(atx, idFusions).execute();
                        }
                    }

                    if (fusions.Bcoord() != null) {
                        for (String bcoord : fusions.Bcoord()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONBCOORD,
                                    MOLECULARMATCHMUTATIONSFUSIONBCOORD.BCOORD,
                                    MOLECULARMATCHMUTATIONSFUSIONBCOORD.MOLECULARMATCHMUTATIONSFUSIONID)
                                    .values(bcoord, idFusions)
                                    .execute();
                        }
                    }

                    if (fusions.Bgreg() != null) {
                        for (MolecularMatchBreg breg : fusions.Bgreg()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONBREG,
                                    MOLECULARMATCHMUTATIONSFUSIONBREG.NUM,
                                    MOLECULARMATCHMUTATIONSFUSIONBREG.TYPE,
                                    MOLECULARMATCHMUTATIONSFUSIONBREG.MOLECULARMATCHMUTATIONSFUSIONID)
                                    .values(breg.num(), breg.type(), idFusions)
                                    .execute();
                        }
                    }

                    if (fusions.Agreg() != null) {
                        for (MolecularMatchAgreg areg : fusions.Agreg()) {
                            context.insertInto(MOLECULARMATCHMUTATIONSFUSIONAREG,
                                    MOLECULARMATCHMUTATIONSFUSIONAREG.NUM,
                                    MOLECULARMATCHMUTATIONSFUSIONAREG.TYPE,
                                    MOLECULARMATCHMUTATIONSFUSIONAREG.MOLECULARMATCHMUTATIONSFUSIONID)
                                    .values(areg.num(), areg.type(), idFusions)
                                    .execute();
                        }
                    }
                }
            }

            MolecularMatchExonsInfo info = mutations.exonsInfo();
            if (info != null) {
                int idExons = context.insertInto(MOLECULARMATCHMUTATIONSEXONSINFO,
                        MOLECULARMATCHMUTATIONSEXONSINFO.TXSTART,
                        MOLECULARMATCHMUTATIONSEXONSINFO.CDSEND,
                        MOLECULARMATCHMUTATIONSEXONSINFO.CHR,
                        MOLECULARMATCHMUTATIONSEXONSINFO.CDSSTART,
                        MOLECULARMATCHMUTATIONSEXONSINFO.TRANSCRIPT,
                        MOLECULARMATCHMUTATIONSEXONSINFO.TXEND,
                        MOLECULARMATCHMUTATIONSEXONSINFO.MOLECULARMATCHMUTATIONSID)
                        .values(info.txStart(), info.cdsEnd(), info.chr(), info.cdsStart(), info.transcript(), info.txEnd(), idMutations)
                        .returning(MOLECULARMATCHMUTATIONSEXONSINFO.ID)
                        .fetchOne()
                        .getValue(MOLECULARMATCHMUTATIONSEXONSINFO.ID);

                int idBoundries = context.insertInto(MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIES,
                        MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIES.MOLECULARMATCHMUTATIONSEXONSINFOID)
                        .values(idExons)
                        .returning(MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIES.ID)
                        .fetchOne()
                        .getValue(MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIES.ID);

                int idExonPosities = context.insertInto(MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIESEXONPOSITIES,
                        MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIESEXONPOSITIES.MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIESID)
                        .values(idBoundries)
                        .returning(MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIESEXONPOSITIES.ID)
                        .fetchOne()
                        .getValue(MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIESEXONPOSITIES.ID);

                if (info.exonBoundaries().exon1() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON1,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON1.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON1.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON1.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon1().start(), info.exonBoundaries().exon1().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon2() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON2,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON2.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON2.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON2.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon2().start(), info.exonBoundaries().exon2().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon3() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON3,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON3.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON3.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON3.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon3().start(), info.exonBoundaries().exon3().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon4() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON4,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON4.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON4.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON4.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon4().start(), info.exonBoundaries().exon4().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon5() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON5,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON5.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON5.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON5.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon5().start(), info.exonBoundaries().exon5().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon6() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON6,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON6.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON6.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON6.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon6().start(), info.exonBoundaries().exon6().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon7() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON7,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON7.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON7.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON7.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon7().start(), info.exonBoundaries().exon7().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon8() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON8,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON8.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON8.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON8.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon8().start(), info.exonBoundaries().exon8().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon9() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON9,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON9.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON9.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON9.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon9().start(), info.exonBoundaries().exon9().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon10() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON10,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON10.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON10.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON10.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon10().start(), info.exonBoundaries().exon10().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon11() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON11,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON11.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON11.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON11.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon11().start(), info.exonBoundaries().exon11().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon12() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON12,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON12.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON12.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON12.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon12().start(), info.exonBoundaries().exon12().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon13() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON13,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON13.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON13.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON13.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon13().start(), info.exonBoundaries().exon13().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon14() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON14,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON14.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON14.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON14.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon14().start(), info.exonBoundaries().exon14().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon15() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON15,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON15.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON15.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON15.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon15().start(), info.exonBoundaries().exon15().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon16() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON16,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON16.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON16.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON16.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon16().start(), info.exonBoundaries().exon16().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon17() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON17,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON17.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON17.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON17.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon17().start(), info.exonBoundaries().exon17().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon18() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON18,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON18.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON18.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON18.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon18().start(), info.exonBoundaries().exon18().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon19() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON19,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON19.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON19.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON19.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon19().start(), info.exonBoundaries().exon19().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon20() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON20,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON20.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON20.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON20.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon20().start(), info.exonBoundaries().exon20().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon21() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON21,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON21.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON21.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON21.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon21().start(), info.exonBoundaries().exon21().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon22() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON22,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON22.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON22.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON22.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon22().start(), info.exonBoundaries().exon22().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon23() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON23,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON23.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON23.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON23.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon23().start(), info.exonBoundaries().exon23().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon24() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON24,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON24.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON24.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON24.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon24().start(), info.exonBoundaries().exon24().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon25() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON25,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON25.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON25.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON25.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon25().start(), info.exonBoundaries().exon25().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon26() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON26,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON26.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON26.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON26.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon26().start(), info.exonBoundaries().exon26().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon27() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON27,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON27.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON27.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON27.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon27().start(), info.exonBoundaries().exon27().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon28() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON28,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON28.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON28.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON28.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon28().start(), info.exonBoundaries().exon28().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon29() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON29,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON29.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON29.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON29.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon29().start(), info.exonBoundaries().exon29().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon30() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON30,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON30.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON30.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON30.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon30().start(), info.exonBoundaries().exon30().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon31() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON31,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON31.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON31.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON31.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon31().start(), info.exonBoundaries().exon31().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon32() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON32,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON32.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON32.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON32.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon32().start(), info.exonBoundaries().exon32().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon33() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON33,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON33.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON33.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON33.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon33().start(), info.exonBoundaries().exon33().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon34() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON34,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON34.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON34.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON34.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon34().start(), info.exonBoundaries().exon34().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon35() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON35,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON35.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON35.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON35.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon35().start(), info.exonBoundaries().exon35().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon36() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON36,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON36.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON36.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON36.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon36().start(), info.exonBoundaries().exon36().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon37() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON37,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON37.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON37.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON37.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon37().start(), info.exonBoundaries().exon37().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon38() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON38,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON38.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON38.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON38.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon38().start(), info.exonBoundaries().exon38().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon39() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON39,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON39.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON39.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON39.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon39().start(), info.exonBoundaries().exon39().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon40() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON40,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON40.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON40.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON40.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon40().start(), info.exonBoundaries().exon40().stop(), idExonPosities)
                            .execute();
                }
                if (info.exonBoundaries().exon41() != null) {
                    context.insertInto(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON41,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON41.START,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON41.END,
                            MOLECULARMATCHMUTATIONSEXONPOSITIESEXON41.MOLECULARMATCHMUTATIONSEXONPOSITIESID)
                            .values(info.exonBoundaries().exon41().start(), info.exonBoundaries().exon41().stop(), idExonPosities)
                            .execute();
                }

            }
        }
    }

    static void deleteAll(@NotNull DSLContext context) {
        // TODO Order from branch to root to avoid constraint violation.
        context.deleteFrom(MOLECULARMATCH).execute();
        context.deleteFrom(MOLECULARMATCHAST).execute();
        context.deleteFrom(MOLECULARMATCHASTLEFT).execute();
        context.deleteFrom(MOLECULARMATCHASTRIGHT).execute();
        context.deleteFrom(MOLECULARMATCHASTRIGHTRIGHT).execute();
        context.deleteFrom(MOLECULARMATCHASTRIGHTLEFT).execute();
        context.deleteFrom(MOLECULARMATCHASTRIGHTLEFTRIGHT).execute();
        context.deleteFrom(MOLECULARMATCHASTRIGHTLEFTLEFT).execute();
        context.deleteFrom(MOLECULARMATCHINSTUTITION).execute();
        context.deleteFrom(MOLECULARMATCHINCLUDEGENE0).execute();
        context.deleteFrom(MOLECULARMATCHEXTERNALID).execute();
        context.deleteFrom(MOLECULARMATCHINCLUDESTAGE0).execute();
        context.deleteFrom(MOLECULARMATCHINCLUDEDRUG1).execute();
        context.deleteFrom(MOLECULARMATCHINCLUDECONDITION1).execute();
        context.deleteFrom(MOLECULARMATCHINCLUDEMUTATION1).execute();
        context.deleteFrom(MOLECULARMATCHINCLUDECONDITION0).execute();
        context.deleteFrom(MOLECULARMATCHINCLUDEMUTATION0).execute();
        context.deleteFrom(MOLECULARMATCHCRITERIAMET).execute();
        context.deleteFrom(MOLECULARMATCHSOURCE).execute();
        context.deleteFrom(MOLECULARMATCHTIEREXPLANATION).execute();
        context.deleteFrom(MOLECULARMATCHTHERAPEUTICCONTEXT).execute();
        context.deleteFrom(MOLECULARMATCHTAGS).execute();
        context.deleteFrom(MOLECULARMATCHVARIANTINFO).execute();
        context.deleteFrom(MOLECULARMATCHVARIANTINFOCONSEQUENCES).execute();
        context.deleteFrom(MOLECULARMATCHVARIANTINFOFUSIONS).execute();
        context.deleteFrom(MOLECULARMATCHVARIANTINFOLOCATIONS).execute();
        context.deleteFrom(MOLECULARMATCHVARIANTINFOLOCATIONSEXONNUMBER).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATION).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONEND).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONSTART).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONCHR).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONPATHOLOGY).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONREF).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONEXON).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONALT).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONSOURCES).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONTRANSCRIPTS).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONCOSMICID).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONDBSNP).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONPOPFREQMAX).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONEXONICFUNC).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONNUCLEOTIDECHANGE).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONPARENTS).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONPARENTSTRANSCRIPTS).execute();
        context.deleteFrom(MOLECULARMATCHCRITERIAUNMET).execute();
        context.deleteFrom(MOLECULARMATCHPREFELANCE).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONS).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSMUTATIONTYPE).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSSOURCE).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSSYNONYMS).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSPATHOLOGY).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSCDNA).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCES).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSTRANSCRIPTCONSEQUENCESEXONNUMBER).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSPARENTS).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSPARENTSTRANSCRIPT).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSWGSDATA).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSWGSDATACLINVARDIS).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSWGSDATACLINVARSIG).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSWGSDATACLINVARSTATUS).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSWGSDATACLINVARDBID).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSWGSDATAFULLAA).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSWGSDATAGENE).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSWGSMAP).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSWGSMAPSYNONYMS).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSWGSMAPPROTCOORDS).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSGRCH37LOCATION).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSGRCH37LOCATIONCONSEQUENCES).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSGRCH37LOCCONSEQUENCESTXSITES).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSGRCH37LOCCONSEQUENCESEXONNUMBER).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSION).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONBCHR).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONAGENE).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONBTX).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONACHR).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONINS).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONBGENE).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONACOORD).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONBORI).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONABAND).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONBBAND).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONAORI).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONATX).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONBCOORD).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONBREG).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSFUSIONAREG).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONSINFO).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIES).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONSINFOBOUNDRIESEXONPOSITIES).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON1).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON2).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON3).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON4).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON5).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON6).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON7).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON8).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON9).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON10).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON11).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON12).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON13).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON14).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON15).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON16).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON17).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON18).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON19).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON20).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON21).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON22).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON23).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON24).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON25).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON26).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON27).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON28).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON29).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON30).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON31).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON32).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON33).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON34).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON35).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON36).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON37).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON38).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON39).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON40).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSEXONPOSITIESEXON41).execute();
    }
}
