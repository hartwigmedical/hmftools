package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCH;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHAST;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHASTLEFT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHASTLEFTLEFT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHASTLEFTRIGHT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHASTRIGHT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHASTRIGHTLEFT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHASTRIGHTRIGHT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONALT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONCHROMOSOME;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONCOSMICID;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONDBSNP;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONEND;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONEXON;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONEXONICFUNC;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONNUCLEOTIDECHANGE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONPARENT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONPARENTTRANSCRIPT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONPATHOLOGY;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONPOPFREQMAX;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONREF;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONSOURCE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONSTART;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCLASSIFICATIONTRANSCRIPT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCRITERIAMET;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHCRITERIAUNMET;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHEXTERNALID;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINCLUDECONDITION0;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINCLUDECONDITION1;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINCLUDEDRUG1;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINCLUDEDRUGCLASS1;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINCLUDEFINDING1;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINCLUDEGENE0;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINCLUDEGENE1;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINCLUDEMUTATION0;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINCLUDEMUTATION1;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINCLUDERESISTANCE1;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINCLUDESTAGE0;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHINSTITUTION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONCDNA;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFO;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON1;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON10;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON11;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON12;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON13;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON14;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON15;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON16;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON17;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON18;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON19;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON2;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON20;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON21;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON22;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON23;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON24;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON25;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON26;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON27;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON28;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON29;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON3;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON30;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON31;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON32;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON33;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON34;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON35;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON36;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON37;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON38;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON39;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON4;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON40;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON41;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON5;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON6;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON7;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON8;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON9;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONFUSION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONFUSIONABAND;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONFUSIONACHROMOSOME;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONFUSIONACOORD;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONFUSIONAGENE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONFUSIONAGENOMICREGION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONFUSIONAORIENTATION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONFUSIONATRANSCRIPT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONFUSIONBBAND;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONFUSIONBCHROMOSOME;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONFUSIONBCOORD;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONFUSIONBGENE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONFUSIONBGENOMICREGION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONFUSIONBORIENTATION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONFUSIONBTRANSCRIPT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONFUSIONINSERT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONGRCH37LOC;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONGRCH37LOCCONSEQUENCE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONGRCH37LOCCONSEQUENCEEXONNUMBER;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONGRCH37LOCCONSEQUENCETXSITE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONMUTATIONTYPE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONPARENT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONPARENTTRANSCRIPT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONPATHOLOGY;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSOURCE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONSYNONYM;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCEEXONNUMBER;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONWGSALOCATION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONWGSALOCATIONCLINVARDBID;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONWGSALOCATIONCLINVARDISEASE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONWGSALOCATIONCLINVARSIG;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONWGSALOCATIONCLINVARSTATUS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONWGSALOCATIONFULLAA;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONWGSALOCATIONGENE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONWGSAMAP;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONWGSAMAPPROTCOORD;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHMUTATIONWGSAMAPSYNONYM;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHPREVALENCE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHSOURCE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTAG;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTHERAPEUTICCONTEXT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTIEREXPLANATION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHVARIANTINFO;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHVARIANTINFOCONSEQUENCE;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHVARIANTINFOFUSION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHVARIANTINFOLOCATION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHVARIANTINFOLOCATIONEXONNUMBER;

import java.util.List;

import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatch;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAst;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAstLeft;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAstLeftLeft;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAstLeftRight;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAstRight;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAstRightLeft;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAstRightRight;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchClassification;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchCriteriaUnmet;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchExonsInfo;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchFusion;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchFusionData;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchFusionGenomicRegion;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchGRCh37Location;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchGRCh37TranscriptConsequence;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchLocation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchMutation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchParent;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchPosition;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchPrevalence;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchSource;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchTag;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchTherapeuticContext;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchTierExplanation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchTranscriptConsequence;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchVariantInfo;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchWGSALocation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchWGSAMap;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

final class MolecularMatchDAOFunctions {

    private MolecularMatchDAOFunctions() {
    }

    static void write(@NotNull DSLContext context, int viccEntryId, @NotNull MolecularMatch molecularMatch) {
        int id = context.insertInto(MOLECULARMATCH,
                MOLECULARMATCH.DIRECTION,
                MOLECULARMATCH.BIOMARKERCLASS,
                MOLECULARMATCH.SCORE,
                MOLECULARMATCH.CLINICALSIGNIFICANCE,
                MOLECULARMATCH.TIER,
                MOLECULARMATCH.AMPCAP,
                MOLECULARMATCH.CIVICVALUE,
                MOLECULARMATCH.REGULATORYBODY,
                MOLECULARMATCH.REGULATORYBODYAPPROVED,
                MOLECULARMATCH.GUIDELINEBODY,
                MOLECULARMATCH.GUIDELINEVERSION,
                MOLECULARMATCH.NOTHERAPYAVAILABLE,
                MOLECULARMATCH.SIXTIER,
                MOLECULARMATCH.MVLD,
                MOLECULARMATCH.AUTOGENERATENARRATIVE,
                MOLECULARMATCH.NARRATIVE,
                MOLECULARMATCH.EXPRESSION,
                MOLECULARMATCH.CUSTOMER,
                MOLECULARMATCH.VERSION,
                MOLECULARMATCH.IDMOLECULARMATCH,
                MOLECULARMATCH.UNIQUEKEY,
                MOLECULARMATCH.HASHKEY,
                MOLECULARMATCH.VICCENTRYID)
                .values(molecularMatch.direction(),
                        molecularMatch.biomarkerClass(),
                        molecularMatch.score(),
                        molecularMatch.clinicalSignificance(),
                        molecularMatch.tier(),
                        molecularMatch.ampcap(),
                        molecularMatch.civicValue(),
                        molecularMatch.regulatoryBody(),
                        molecularMatch.regulatoryBodyApproved(),
                        molecularMatch.guidelineBody(),
                        molecularMatch.guidelineVersion(),
                        molecularMatch.noTherapyAvailable(),
                        molecularMatch.sixtier(),
                        molecularMatch.mvld(),
                        molecularMatch.autoGenerateNarrative(),
                        molecularMatch.narrative(),
                        molecularMatch.expression(),
                        molecularMatch.customer(),
                        molecularMatch.version(),
                        molecularMatch.id(),
                        molecularMatch.uniqueKey(),
                        molecularMatch.hashKey(),
                        viccEntryId)
                .returning(MOLECULARMATCH.ID)
                .fetchOne()
                .getValue(MOLECULARMATCH.ID);

        insertMutations(context, molecularMatch.mutations(), id);
        insertVariantInfos(context, molecularMatch.variantInfos(), id);
        insertAst(context, molecularMatch.ast(), id);
        insertClassifications(context, molecularMatch.classifications(), id);

        for (String includeGene1 : molecularMatch.includeGene1()) {
            context.insertInto(MOLECULARMATCHINCLUDEGENE1,
                    MOLECULARMATCHINCLUDEGENE1.INCLUDEGENE1,
                    MOLECULARMATCHINCLUDEGENE1.MOLECULARMATCHID).values(includeGene1, id).execute();
        }

        for (String includeFinding1 : molecularMatch.includeFinding1()) {
            context.insertInto(MOLECULARMATCHINCLUDEFINDING1,
                    MOLECULARMATCHINCLUDEFINDING1.INCLUDEFINDING1,
                    MOLECULARMATCHINCLUDEFINDING1.MOLECULARMATCHID).values(includeFinding1, id).execute();
        }

        for (String includeCondition1 : molecularMatch.includeCondition1()) {
            context.insertInto(MOLECULARMATCHINCLUDECONDITION1,
                    MOLECULARMATCHINCLUDECONDITION1.INCLUDECONDITION1,
                    MOLECULARMATCHINCLUDECONDITION1.MOLECULARMATCHID).values(includeCondition1, id).execute();
        }

        for (String includeMutation1 : molecularMatch.includeMutation1()) {
            context.insertInto(MOLECULARMATCHINCLUDEMUTATION1,
                    MOLECULARMATCHINCLUDEMUTATION1.INCLUDEMUTATION1,
                    MOLECULARMATCHINCLUDEMUTATION1.MOLECULARMATCHID).values(includeMutation1, id).execute();
        }

        for (String includeDrug1 : molecularMatch.includeDrug1()) {
            context.insertInto(MOLECULARMATCHINCLUDEDRUG1,
                    MOLECULARMATCHINCLUDEDRUG1.INCLUDEDRUG1,
                    MOLECULARMATCHINCLUDEDRUG1.MOLECULARMATCHID).values(includeDrug1, id).execute();
        }

        for (String includeDrugClass1 : molecularMatch.includeDrugClass1()) {
            context.insertInto(MOLECULARMATCHINCLUDEDRUGCLASS1,
                    MOLECULARMATCHINCLUDEDRUGCLASS1.INCLUDEDRUGCLASS1,
                    MOLECULARMATCHINCLUDEDRUGCLASS1.MOLECULARMATCHID).values(includeDrugClass1, id).execute();
        }

        for (String includeResistance1 : molecularMatch.includeResistance1()) {
            context.insertInto(MOLECULARMATCHINCLUDERESISTANCE1,
                    MOLECULARMATCHINCLUDERESISTANCE1.INCLUDERESISTANCE1,
                    MOLECULARMATCHINCLUDERESISTANCE1.MOLECULARMATCHID).values(includeResistance1, id).execute();
        }

        for (String includeStage0 : molecularMatch.includeStage0()) {
            context.insertInto(MOLECULARMATCHINCLUDESTAGE0,
                    MOLECULARMATCHINCLUDESTAGE0.INCLUDESTAGE0,
                    MOLECULARMATCHINCLUDESTAGE0.MOLECULARMATCHID).values(includeStage0, id).execute();
        }

        for (String includeGene0 : molecularMatch.includeGene0()) {
            context.insertInto(MOLECULARMATCHINCLUDEGENE0,
                    MOLECULARMATCHINCLUDEGENE0.INCLUDEGENE0,
                    MOLECULARMATCHINCLUDEGENE0.MOLECULARMATCHID).values(includeGene0, id).execute();
        }

        for (String includeCondition0 : molecularMatch.includeCondition0()) {
            context.insertInto(MOLECULARMATCHINCLUDECONDITION0,
                    MOLECULARMATCHINCLUDECONDITION0.INCLUDECONDITION0,
                    MOLECULARMATCHINCLUDECONDITION0.MOLECULARMATCHID).values(includeCondition0, id).execute();
        }

        for (String includeMutation0 : molecularMatch.includeMutation0()) {
            context.insertInto(MOLECULARMATCHINCLUDEMUTATION0,
                    MOLECULARMATCHINCLUDEMUTATION0.INCLUDEMUTATION0,
                    MOLECULARMATCHINCLUDEMUTATION0.MOLECULARMATCHID).values(includeMutation0, id).execute();
        }

        for (String criteriaMet : molecularMatch.criteriaMets()) {
            context.insertInto(MOLECULARMATCHCRITERIAMET, MOLECULARMATCHCRITERIAMET.CRITERIAMET, MOLECULARMATCHCRITERIAMET.MOLECULARMATCHID)
                    .values(criteriaMet, id)
                    .execute();
        }

        for (String institution : molecularMatch.institutions()) {
            context.insertInto(MOLECULARMATCHINSTITUTION, MOLECULARMATCHINSTITUTION.INSTITUTION, MOLECULARMATCHINSTITUTION.MOLECULARMATCHID)
                    .values(institution, id)
                    .execute();
        }

        for (String externalId : molecularMatch.externalIds()) {
            context.insertInto(MOLECULARMATCHEXTERNALID, MOLECULARMATCHEXTERNALID.EXTERNALID, MOLECULARMATCHEXTERNALID.MOLECULARMATCHID)
                    .values(externalId, id)
                    .execute();
        }

        for (MolecularMatchSource source : molecularMatch.sources()) {
            context.insertInto(MOLECULARMATCHSOURCE,
                    MOLECULARMATCHSOURCE.NAME,
                    MOLECULARMATCHSOURCE.TYPE,
                    MOLECULARMATCHSOURCE.SUBTYPE,
                    MOLECULARMATCHSOURCE.VALID,
                    MOLECULARMATCHSOURCE.PUBID,
                    MOLECULARMATCHSOURCE.LINK,
                    MOLECULARMATCHSOURCE.TRIALID,
                    MOLECULARMATCHSOURCE.TRIALPHASE,
                    MOLECULARMATCHSOURCE.YEAR,
                    MOLECULARMATCHSOURCE.FUNCTIONALCONSEQUENCE,
                    MOLECULARMATCHSOURCE.INSTITUTION,
                    MOLECULARMATCHSOURCE.TRUSTRATING,
                    MOLECULARMATCHSOURCE.SUPPRESS,
                    MOLECULARMATCHSOURCE.IDSOURCE,
                    MOLECULARMATCHSOURCE.MOLECULARMATCHID)
                    .values(source.name(),
                            source.type(),
                            source.subType(),
                            source.valid(),
                            source.pubId(),
                            source.link(),
                            source.trialId(),
                            source.trialPhase(),
                            source.year(),
                            source.functionalConsequence(),
                            source.institution(),
                            source.trustRating(),
                            source.suppress(),
                            source.id(),
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
                    MOLECULARMATCHTHERAPEUTICCONTEXT.NAME,
                    MOLECULARMATCHTHERAPEUTICCONTEXT.FACET,
                    MOLECULARMATCHTHERAPEUTICCONTEXT.SUPPRESS,
                    MOLECULARMATCHTHERAPEUTICCONTEXT.VALID,
                    MOLECULARMATCHTHERAPEUTICCONTEXT.MOLECULARMATCHID)
                    .values(therapeuticContext.name(),
                            therapeuticContext.facet(),
                            therapeuticContext.suppress(),
                            therapeuticContext.valid(),
                            id)
                    .execute();
        }

        for (MolecularMatchTag tag : molecularMatch.tags()) {
            context.insertInto(MOLECULARMATCHTAG,
                    MOLECULARMATCHTAG.TERM,
                    MOLECULARMATCHTAG.FACET,
                    MOLECULARMATCHTAG.FILTERTYPE,
                    MOLECULARMATCHTAG.PRIORITY,
                    MOLECULARMATCHTAG.TRANSCRIPT,
                    MOLECULARMATCHTAG.VALID,
                    MOLECULARMATCHTAG.GENERATEDBY,
                    MOLECULARMATCHTAG.GENERATEDBYTERM,
                    MOLECULARMATCHTAG.ISNEW,
                    MOLECULARMATCHTAG.PRIMARYVALUE,
                    MOLECULARMATCHTAG.CUSTOM,
                    MOLECULARMATCHTAG.SUPPRESS,
                    MOLECULARMATCHTAG.MANUALSUPPRESS,
                    MOLECULARMATCHTAG.COMPOSITE,
                    MOLECULARMATCHTAG.COMPOSITEKEY,
                    MOLECULARMATCHTAG.MOLECULARMATCHID)
                    .values(tag.term(),
                            tag.facet(),
                            tag.filterType(),
                            tag.priority(),
                            tag.transcript(),
                            tag.valid(),
                            tag.generatedBy(),
                            tag.generatedByTerm(),
                            tag.isNew(),
                            tag.primary(),
                            tag.custom(),
                            tag.suppress(),
                            tag.manualSuppress(),
                            tag.composite(),
                            tag.compositeKey(),
                            id)
                    .execute();
        }

        for (MolecularMatchCriteriaUnmet criteriaUnmet : molecularMatch.criteriaUnmets()) {
            context.insertInto(MOLECULARMATCHCRITERIAUNMET,
                    MOLECULARMATCHCRITERIAUNMET.TERM,
                    MOLECULARMATCHCRITERIAUNMET.FILTERTYPE,
                    MOLECULARMATCHCRITERIAUNMET.PRIORITY,
                    MOLECULARMATCHCRITERIAUNMET.FACET,
                    MOLECULARMATCHCRITERIAUNMET.VALID,
                    MOLECULARMATCHCRITERIAUNMET.TRANSCRIPT,
                    MOLECULARMATCHCRITERIAUNMET.ISNEW,
                    MOLECULARMATCHCRITERIAUNMET.GENERATEDBY,
                    MOLECULARMATCHCRITERIAUNMET.GENERATEDBYTERM,
                    MOLECULARMATCHCRITERIAUNMET.SUPPRESS,
                    MOLECULARMATCHCRITERIAUNMET.MANUALSUPPRESS,
                    MOLECULARMATCHCRITERIAUNMET.PRIMARYVALUE,
                    MOLECULARMATCHCRITERIAUNMET.COMPOSITEKEY,
                    MOLECULARMATCHCRITERIAUNMET.CUSTOM,
                    MOLECULARMATCHCRITERIAUNMET.MOLECULARMATCHID)
                    .values(criteriaUnmet.term(),
                            criteriaUnmet.filterType(),
                            criteriaUnmet.priority(),
                            criteriaUnmet.facet(),
                            criteriaUnmet.valid(),
                            criteriaUnmet.transcript(),
                            criteriaUnmet.isNew(),
                            criteriaUnmet.generatedBy(),
                            criteriaUnmet.generatedByTerm(),
                            criteriaUnmet.suppress(),
                            criteriaUnmet.manualSuppress(),
                            criteriaUnmet.primary(),
                            criteriaUnmet.compositeKey(),
                            criteriaUnmet.custom(),
                            id)
                    .execute();
        }

        for (MolecularMatchPrevalence prevalence : molecularMatch.prevalences()) {
            context.insertInto(MOLECULARMATCHPREVALENCE,
                    MOLECULARMATCHPREVALENCE.STUDYID,
                    MOLECULARMATCHPREVALENCE.COUNT,
                    MOLECULARMATCHPREVALENCE.SAMPLES,
                    MOLECULARMATCHPREVALENCE.PERCENT,
                    MOLECULARMATCHPREVALENCE.MOLECULAR,
                    MOLECULARMATCHPREVALENCE.CONDITIONVALUE,
                    MOLECULARMATCHPREVALENCE.MOLECULARMATCHID)
                    .values(prevalence.studyId(),
                            prevalence.count(),
                            prevalence.samples(),
                            prevalence.percent(),
                            prevalence.molecular(),
                            prevalence.condition(),
                            id)
                    .execute();
        }
    }

    private static void insertAst(@NotNull DSLContext context, @NotNull MolecularMatchAst ast, int molecularMatchId) {
        int astId = context.insertInto(MOLECULARMATCHAST,
                MOLECULARMATCHAST.TYPE,
                MOLECULARMATCHAST.RAW,
                MOLECULARMATCHAST.VALUE,
                MOLECULARMATCHAST.OPERATOR,
                MOLECULARMATCHAST.MOLECULARMATCHID)
                .values(ast.type(), ast.raw(), ast.value(), ast.operator(), molecularMatchId)
                .returning(MOLECULARMATCHAST.ID)
                .fetchOne()
                .getValue(MOLECULARMATCHAST.ID);

        MolecularMatchAstLeft astLeft = ast.left();
        if (astLeft != null) {
            int astLeftId = context.insertInto(MOLECULARMATCHASTLEFT,
                    MOLECULARMATCHASTLEFT.TYPE,
                    MOLECULARMATCHASTLEFT.RAW,
                    MOLECULARMATCHASTLEFT.VALUE,
                    MOLECULARMATCHASTLEFT.OPERATOR,
                    MOLECULARMATCHASTLEFT.MOLECULARMATCHASTID)
                    .values(astLeft.type(), astLeft.raw(), astLeft.value(), astLeft.operator(), astId)
                    .returning(MOLECULARMATCHASTLEFT.ID)
                    .fetchOne()
                    .getValue(MOLECULARMATCHASTLEFT.ID);

            MolecularMatchAstLeftLeft astLeftLeft = astLeft.left();
            if (astLeftLeft != null) {
                context.insertInto(MOLECULARMATCHASTLEFTLEFT,
                        MOLECULARMATCHASTLEFTLEFT.TYPE,
                        MOLECULARMATCHASTLEFTLEFT.RAW,
                        MOLECULARMATCHASTLEFTLEFT.VALUE,
                        MOLECULARMATCHASTLEFTLEFT.OPERATOR,
                        MOLECULARMATCHASTLEFTLEFT.MOLECULARMATCHASTLEFTID)
                        .values(astLeftLeft.type(), astLeftLeft.raw(), astLeftLeft.value(), astLeftLeft.operator(), astLeftId)
                        .execute();
            }

            MolecularMatchAstLeftRight astLeftRight = astLeft.right();
            if (astLeftRight != null) {
                context.insertInto(MOLECULARMATCHASTLEFTRIGHT,
                        MOLECULARMATCHASTLEFTRIGHT.TYPE,
                        MOLECULARMATCHASTLEFTRIGHT.RAW,
                        MOLECULARMATCHASTLEFTRIGHT.VALUE,
                        MOLECULARMATCHASTLEFTRIGHT.OPERATOR,
                        MOLECULARMATCHASTLEFTRIGHT.MOLECULARMATCHASTLEFTID)
                        .values(astLeftRight.type(), astLeftRight.raw(), astLeftRight.value(), astLeftRight.operator(), astLeftId)
                        .execute();
            }
        }

        MolecularMatchAstRight astRight = ast.right();
        if (astRight != null) {
            int astRightId = context.insertInto(MOLECULARMATCHASTRIGHT,
                    MOLECULARMATCHASTRIGHT.TYPE,
                    MOLECULARMATCHASTRIGHT.RAW,
                    MOLECULARMATCHASTRIGHT.VALUE,
                    MOLECULARMATCHASTRIGHT.OPERATOR,
                    MOLECULARMATCHASTRIGHT.MOLECULARMATCHASTID)
                    .values(astRight.type(), astRight.raw(), astRight.value(), astRight.operator(), astId)
                    .returning(MOLECULARMATCHASTRIGHT.ID)
                    .fetchOne()
                    .getValue(MOLECULARMATCHASTRIGHT.ID);

            MolecularMatchAstRightLeft astRightLeft = astRight.left();
            if (astRightLeft != null) {
                context.insertInto(MOLECULARMATCHASTRIGHTLEFT,
                        MOLECULARMATCHASTRIGHTLEFT.TYPE,
                        MOLECULARMATCHASTRIGHTLEFT.RAW,
                        MOLECULARMATCHASTRIGHTLEFT.VALUE,
                        MOLECULARMATCHASTRIGHTLEFT.OPERATOR,
                        MOLECULARMATCHASTRIGHTLEFT.MOLECULARMATCHASTRIGHTID)
                        .values(astRightLeft.type(), astRightLeft.raw(), astRightLeft.value(), astRightLeft.operator(), astRightId)
                        .execute();
            }

            MolecularMatchAstRightRight astRightRight = astRight.right();
            if (astRightRight != null) {
                context.insertInto(MOLECULARMATCHASTRIGHTRIGHT,
                        MOLECULARMATCHASTRIGHTRIGHT.TYPE,
                        MOLECULARMATCHASTRIGHTRIGHT.RAW,
                        MOLECULARMATCHASTRIGHTRIGHT.VALUE,
                        MOLECULARMATCHASTRIGHTRIGHT.OPERATOR,
                        MOLECULARMATCHASTRIGHTRIGHT.MOLECULARMATCHASTRIGHTID)
                        .values(astRightRight.type(), astRightRight.raw(), astRightRight.value(), astRightRight.operator(), astRightId)
                        .execute();
            }

        }
    }

    private static void insertVariantInfos(@NotNull DSLContext context, @NotNull List<MolecularMatchVariantInfo> variantInfos,
            int molecularMatchId) {
        for (MolecularMatchVariantInfo variantInfo : variantInfos) {
            int variantInfoId = context.insertInto(MOLECULARMATCHVARIANTINFO,
                    MOLECULARMATCHVARIANTINFO.NAME,
                    MOLECULARMATCHVARIANTINFO.GENE,
                    MOLECULARMATCHVARIANTINFO.TRANSCRIPT,
                    MOLECULARMATCHVARIANTINFO.CLASSIFICATION,
                    MOLECULARMATCHVARIANTINFO.GENEFUSIONPARTNER,
                    MOLECULARMATCHVARIANTINFO.COSMICID,
                    MOLECULARMATCHVARIANTINFO.POPFREQMAX,
                    MOLECULARMATCHVARIANTINFO.MOLECULARMATCHID)
                    .values(variantInfo.name(),
                            variantInfo.gene(),
                            variantInfo.transcript(),
                            variantInfo.classification(),
                            variantInfo.geneFusionPartner(),
                            variantInfo.cosmicId(),
                            variantInfo.popFreqMax(),
                            molecularMatchId)
                    .returning(MOLECULARMATCHVARIANTINFO.ID)
                    .fetchOne()
                    .getValue(MOLECULARMATCHVARIANTINFO.ID);

            for (String consequence : variantInfo.consequences()) {
                context.insertInto(MOLECULARMATCHVARIANTINFOCONSEQUENCE,
                        MOLECULARMATCHVARIANTINFOCONSEQUENCE.CONSEQUENCE,
                        MOLECULARMATCHVARIANTINFOCONSEQUENCE.MOLECULARMATCHVARIANTINFOID).values(consequence, variantInfoId).execute();
            }

            for (MolecularMatchFusion fusion : variantInfo.fusions()) {
                context.insertInto(MOLECULARMATCHVARIANTINFOFUSION,
                        MOLECULARMATCHVARIANTINFOFUSION.CHR,
                        MOLECULARMATCHVARIANTINFOFUSION.REFERENCEGENOME,
                        MOLECULARMATCHVARIANTINFOFUSION.LBPWREP,
                        MOLECULARMATCHVARIANTINFOFUSION.LBPWLEP,
                        MOLECULARMATCHVARIANTINFOFUSION.RBPWREP,
                        MOLECULARMATCHVARIANTINFOFUSION.RBPWLEP,
                        MOLECULARMATCHVARIANTINFOFUSION.INTRONNUMBER,
                        MOLECULARMATCHVARIANTINFOFUSION.EXONNUMBER,
                        MOLECULARMATCHVARIANTINFOFUSION.MOLECULARMATCHVARIANTINFOID)
                        .values(fusion.chr(),
                                fusion.referenceGenome(),
                                fusion.LBPWREP(),
                                fusion.LBPWLEP(),
                                fusion.RBPWREP(),
                                fusion.RBPWLEP(),
                                fusion.intronNumber(),
                                fusion.exonNumber(),
                                variantInfoId)
                        .execute();
            }

            for (MolecularMatchLocation location : variantInfo.locations()) {
                int locationId = context.insertInto(MOLECULARMATCHVARIANTINFOLOCATION,
                        MOLECULARMATCHVARIANTINFOLOCATION.CHR,
                        MOLECULARMATCHVARIANTINFOLOCATION.START,
                        MOLECULARMATCHVARIANTINFOLOCATION.STOP,
                        MOLECULARMATCHVARIANTINFOLOCATION.REF,
                        MOLECULARMATCHVARIANTINFOLOCATION.ALT,
                        MOLECULARMATCHVARIANTINFOLOCATION.CDNA,
                        MOLECULARMATCHVARIANTINFOLOCATION.AMINOACIDCHANGE,
                        MOLECULARMATCHVARIANTINFOLOCATION.REFERENCEGENOME,
                        MOLECULARMATCHVARIANTINFOLOCATION.STRAND,
                        MOLECULARMATCHVARIANTINFOLOCATION.INTRONNUMBER,
                        MOLECULARMATCHVARIANTINFOLOCATION.MOLECULARMATCHVARIANTINFOID)
                        .values(location.chr(),
                                location.start(),
                                location.stop(),
                                location.ref(),
                                location.alt(),
                                location.cdna(),
                                location.aminoAcidChange(),
                                location.referenceGenome(),
                                location.strand(),
                                location.intronNumber(),
                                variantInfoId)
                        .returning(MOLECULARMATCHVARIANTINFOLOCATION.ID)
                        .fetchOne()
                        .getValue(MOLECULARMATCHVARIANTINFOLOCATION.ID);

                for (String exonNumber : location.exonNumbers()) {
                    context.insertInto(MOLECULARMATCHVARIANTINFOLOCATIONEXONNUMBER,
                            MOLECULARMATCHVARIANTINFOLOCATIONEXONNUMBER.EXONNUMBER,
                            MOLECULARMATCHVARIANTINFOLOCATIONEXONNUMBER.MOLECULARMATCHVARIANTINFOLOCATIONID)
                            .values(exonNumber, locationId)
                            .execute();
                }
            }
        }
    }

    private static void insertClassifications(@NotNull DSLContext context, @NotNull List<MolecularMatchClassification> classifications,
            int molecularMatchId) {
        for (MolecularMatchClassification classification : classifications) {
            int classificationId = context.insertInto(MOLECULARMATCHCLASSIFICATION,
                    MOLECULARMATCHCLASSIFICATION.NAME,
                    MOLECULARMATCHCLASSIFICATION.GENESYMBOL,
                    MOLECULARMATCHCLASSIFICATION.EXPANDGENESEARCH,
                    MOLECULARMATCHCLASSIFICATION.TRANSCRIPT,
                    MOLECULARMATCHCLASSIFICATION.CLASSIFICATION,
                    MOLECULARMATCHCLASSIFICATION.CLASSIFICATIONOVERRIDE,
                    MOLECULARMATCHCLASSIFICATION.COPYNUMBERTYPE,
                    MOLECULARMATCHCLASSIFICATION.DRUGSAPPROVEDONLABELCOUNT,
                    MOLECULARMATCHCLASSIFICATION.DRUGSAPPROVEDOFFLABELCOUNT,
                    MOLECULARMATCHCLASSIFICATION.DRUGSEXPERIMENTALCOUNT,
                    MOLECULARMATCHCLASSIFICATION.TRIALCOUNT,
                    MOLECULARMATCHCLASSIFICATION.PUBLICATIONCOUNT,
                    MOLECULARMATCHCLASSIFICATION.ROOTTERM,
                    MOLECULARMATCHCLASSIFICATION.ALIAS,
                    MOLECULARMATCHCLASSIFICATION.PRIORITY,
                    MOLECULARMATCHCLASSIFICATION.DESCRIPTION,
                    MOLECULARMATCHCLASSIFICATION.MOLECULARMATCHID)
                    .values(classification.name(),
                            classification.geneSymbol(),
                            classification.expandGeneSearch(),
                            classification.transcript(),
                            classification.classification(),
                            classification.classificationOverride(),
                            classification.copyNumberType(),
                            classification.drugsApprovedOnLabelCount(),
                            classification.drugsApprovedOffLabelCount(),
                            classification.drugsExperimentalCount(),
                            classification.trialCount(),
                            classification.publicationCount(),
                            classification.rootTerm(),
                            classification.alias(),
                            classification.priority(),
                            classification.description(),
                            molecularMatchId)
                    .returning(MOLECULARMATCHCLASSIFICATION.ID)
                    .fetchOne()
                    .getValue(MOLECULARMATCHCLASSIFICATION.ID);

            for (String transcript : classification.transcripts()) {
                context.insertInto(MOLECULARMATCHCLASSIFICATIONTRANSCRIPT,
                        MOLECULARMATCHCLASSIFICATIONTRANSCRIPT.TRANSCRIPT,
                        MOLECULARMATCHCLASSIFICATIONTRANSCRIPT.MOLECULARMATCHCLASSIFICATIONID)
                        .values(transcript, classificationId)
                        .execute();
            }

            for (String chromosome : classification.chromosomes()) {
                context.insertInto(MOLECULARMATCHCLASSIFICATIONCHROMOSOME,
                        MOLECULARMATCHCLASSIFICATIONCHROMOSOME.CHROMOSOME,
                        MOLECULARMATCHCLASSIFICATIONCHROMOSOME.MOLECULARMATCHCLASSIFICATIONID)
                        .values(chromosome, classificationId)
                        .execute();
            }

            for (String start : classification.starts()) {
                context.insertInto(MOLECULARMATCHCLASSIFICATIONSTART,
                        MOLECULARMATCHCLASSIFICATIONSTART.START,
                        MOLECULARMATCHCLASSIFICATIONSTART.MOLECULARMATCHCLASSIFICATIONID).values(start, classificationId).execute();
            }

            for (String end : classification.ends()) {
                context.insertInto(MOLECULARMATCHCLASSIFICATIONEND,
                        MOLECULARMATCHCLASSIFICATIONEND.END,
                        MOLECULARMATCHCLASSIFICATIONEND.MOLECULARMATCHCLASSIFICATIONID).values(end, classificationId).execute();
            }

            for (String ref : classification.refs()) {
                context.insertInto(MOLECULARMATCHCLASSIFICATIONREF,
                        MOLECULARMATCHCLASSIFICATIONREF.REF,
                        MOLECULARMATCHCLASSIFICATIONREF.MOLECULARMATCHCLASSIFICATIONID).values(ref, classificationId).execute();
            }

            for (String alt : classification.alts()) {
                context.insertInto(MOLECULARMATCHCLASSIFICATIONALT,
                        MOLECULARMATCHCLASSIFICATIONALT.ALT,
                        MOLECULARMATCHCLASSIFICATIONALT.MOLECULARMATCHCLASSIFICATIONID).values(alt, classificationId).execute();
            }

            for (String nucleotideChange : classification.nucleotideChanges()) {
                context.insertInto(MOLECULARMATCHCLASSIFICATIONNUCLEOTIDECHANGE,
                        MOLECULARMATCHCLASSIFICATIONNUCLEOTIDECHANGE.NUCLEOTIDECHANGE,
                        MOLECULARMATCHCLASSIFICATIONNUCLEOTIDECHANGE.MOLECULARMATCHCLASSIFICATIONID)
                        .values(nucleotideChange, classificationId)
                        .execute();
            }

            for (String exon : classification.exons()) {
                context.insertInto(MOLECULARMATCHCLASSIFICATIONEXON,
                        MOLECULARMATCHCLASSIFICATIONEXON.EXON,
                        MOLECULARMATCHCLASSIFICATIONEXON.MOLECULARMATCHCLASSIFICATIONID).values(exon, classificationId).execute();
            }

            for (String exonicFunc : classification.exonicFuncs()) {
                context.insertInto(MOLECULARMATCHCLASSIFICATIONEXONICFUNC,
                        MOLECULARMATCHCLASSIFICATIONEXONICFUNC.EXONICFUNC,
                        MOLECULARMATCHCLASSIFICATIONEXONICFUNC.MOLECULARMATCHCLASSIFICATIONID)
                        .values(exonicFunc, classificationId)
                        .execute();
            }

            for (String pathology : classification.pathology()) {
                context.insertInto(MOLECULARMATCHCLASSIFICATIONPATHOLOGY,
                        MOLECULARMATCHCLASSIFICATIONPATHOLOGY.PATHOLOGY,
                        MOLECULARMATCHCLASSIFICATIONPATHOLOGY.MOLECULARMATCHCLASSIFICATIONID).values(pathology, classificationId).execute();
            }

            for (String source : classification.sources()) {
                context.insertInto(MOLECULARMATCHCLASSIFICATIONSOURCE,
                        MOLECULARMATCHCLASSIFICATIONSOURCE.SOURCE,
                        MOLECULARMATCHCLASSIFICATIONSOURCE.MOLECULARMATCHCLASSIFICATIONID).values(source, classificationId).execute();
            }

            for (String dbSNP : classification.dbSNPs()) {
                context.insertInto(MOLECULARMATCHCLASSIFICATIONDBSNP,
                        MOLECULARMATCHCLASSIFICATIONDBSNP.DBSNP,
                        MOLECULARMATCHCLASSIFICATIONDBSNP.MOLECULARMATCHCLASSIFICATIONID).values(dbSNP, classificationId).execute();
            }

            for (String cosmicId : classification.cosmicIds()) {
                context.insertInto(MOLECULARMATCHCLASSIFICATIONCOSMICID,
                        MOLECULARMATCHCLASSIFICATIONCOSMICID.COSMICID,
                        MOLECULARMATCHCLASSIFICATIONCOSMICID.MOLECULARMATCHCLASSIFICATIONID).values(cosmicId, classificationId).execute();
            }

            for (String popFreqMax : classification.popFreqMaxes()) {
                context.insertInto(MOLECULARMATCHCLASSIFICATIONPOPFREQMAX,
                        MOLECULARMATCHCLASSIFICATIONPOPFREQMAX.POPFREQMAX,
                        MOLECULARMATCHCLASSIFICATIONPOPFREQMAX.MOLECULARMATCHCLASSIFICATIONID)
                        .values(popFreqMax, classificationId)
                        .execute();
            }

            for (MolecularMatchParent parent : classification.parents()) {
                int parentId = context.insertInto(MOLECULARMATCHCLASSIFICATIONPARENT,
                        MOLECULARMATCHCLASSIFICATIONPARENT.NAME,
                        MOLECULARMATCHCLASSIFICATIONPARENT.TYPE,
                        MOLECULARMATCHCLASSIFICATIONPARENT.ACTIONABLEPARENT,
                        MOLECULARMATCHCLASSIFICATIONPARENT.MOLECULARMATCHCLASSIFICATIONID)
                        .values(parent.name(), parent.type(), parent.actionableParent(), classificationId)
                        .returning(MOLECULARMATCHCLASSIFICATIONPARENT.ID)
                        .fetchOne()
                        .getValue(MOLECULARMATCHCLASSIFICATIONPARENT.ID);

                for (String transcript : parent.transcripts()) {
                    context.insertInto(MOLECULARMATCHCLASSIFICATIONPARENTTRANSCRIPT,
                            MOLECULARMATCHCLASSIFICATIONPARENTTRANSCRIPT.TRANSCRIPT,
                            MOLECULARMATCHCLASSIFICATIONPARENTTRANSCRIPT.MOLECULARMATCHCLASSIFICATIONPARENTID)
                            .values(transcript, parentId)
                            .execute();
                }
            }
        }
    }

    private static void insertMutations(@NotNull DSLContext context, @NotNull List<MolecularMatchMutation> mutations,
            int molecularMatchId) {
        for (MolecularMatchMutation mutation : mutations) {
            int mutationId = context.insertInto(MOLECULARMATCHMUTATION,
                    MOLECULARMATCHMUTATION.GENESYMBOL,
                    MOLECULARMATCHMUTATION.NAME,
                    MOLECULARMATCHMUTATION.TRANSCRIPTRECOGNIZED,
                    MOLECULARMATCHMUTATION.TRANSCRIPT,
                    MOLECULARMATCHMUTATION.LONGESTTRANSCRIPT,
                    MOLECULARMATCHMUTATION.UNIPROTTRANSCRIPT,
                    MOLECULARMATCHMUTATION.DESCRIPTION,
                    MOLECULARMATCHMUTATION.SRC,
                    MOLECULARMATCHMUTATION.IDMUTATION,
                    MOLECULARMATCHMUTATION.MOLECULARMATCHID)
                    .values(mutation.geneSymbol(),
                            mutation.name(),
                            mutation.transcriptRecognized(),
                            mutation.transcript(),
                            mutation.longestTranscript(),
                            mutation.uniprotTranscript(),
                            mutation.description(),
                            mutation.src(),
                            mutation.id(),
                            molecularMatchId)
                    .returning(MOLECULARMATCHMUTATION.ID)
                    .fetchOne()
                    .getValue(MOLECULARMATCHMUTATION.ID);

            for (String mutationType : mutation.mutationTypes()) {
                context.insertInto(MOLECULARMATCHMUTATIONMUTATIONTYPE,
                        MOLECULARMATCHMUTATIONMUTATIONTYPE.MUTATIONTYPE,
                        MOLECULARMATCHMUTATIONMUTATIONTYPE.MOLECULARMATCHMUTATIONID).values(mutationType, mutationId).execute();
            }

            for (String source : mutation.sources()) {
                context.insertInto(MOLECULARMATCHMUTATIONSOURCE,
                        MOLECULARMATCHMUTATIONSOURCE.SOURCE,
                        MOLECULARMATCHMUTATIONSOURCE.MOLECULARMATCHMUTATIONID).values(source, mutationId).execute();
            }

            for (String synonym : mutation.synonyms()) {
                context.insertInto(MOLECULARMATCHMUTATIONSYNONYM,
                        MOLECULARMATCHMUTATIONSYNONYM.SYNONYM,
                        MOLECULARMATCHMUTATIONSYNONYM.MOLECULARMATCHMUTATIONID).values(synonym, mutationId).execute();
            }

            for (String pathology : mutation.pathology()) {
                context.insertInto(MOLECULARMATCHMUTATIONPATHOLOGY,
                        MOLECULARMATCHMUTATIONPATHOLOGY.PATHOLOGY,
                        MOLECULARMATCHMUTATIONPATHOLOGY.MOLECULARMATCHMUTATIONID).values(pathology, mutationId).execute();
            }

            for (String cDNA : mutation.cDNA()) {
                context.insertInto(MOLECULARMATCHMUTATIONCDNA,
                        MOLECULARMATCHMUTATIONCDNA.CDNA,
                        MOLECULARMATCHMUTATIONCDNA.MOLECULARMATCHMUTATIONID).values(cDNA, mutationId).execute();
            }

            for (MolecularMatchTranscriptConsequence transcriptConsequence : mutation.transcriptConsequences()) {
                int transcriptConsequenceId = context.insertInto(MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCE,
                        MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCE.CHR,
                        MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCE.START,
                        MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCE.STOP,
                        MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCE.REF,
                        MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCE.ALT,
                        MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCE.REFERENCEGENOME,
                        MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCE.TRANSCRIPT,
                        MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCE.STRAND,
                        MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCE.CDNA,
                        MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCE.AMINOACIDCHANGE,
                        MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCE.INTRONNUMBER,
                        MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCE.SUPPRESS,
                        MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCE.CUSTOM,
                        MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCE.VALIDATED,
                        MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCE.COMPOSITEKEY,
                        MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCE.MOLECULARMATCHMUTATIONID)
                        .values(transcriptConsequence.chr(),
                                transcriptConsequence.start(),
                                transcriptConsequence.stop(),
                                transcriptConsequence.ref(),
                                transcriptConsequence.alt(),
                                transcriptConsequence.referenceGenome(),
                                transcriptConsequence.transcript(),
                                transcriptConsequence.strand(),
                                transcriptConsequence.cdna(),
                                transcriptConsequence.aminoAcidChange(),
                                transcriptConsequence.intronNumber(),
                                transcriptConsequence.suppress(),
                                transcriptConsequence.custom(),
                                transcriptConsequence.validated(),
                                transcriptConsequence.compositeKey(),
                                mutationId)
                        .returning(MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCE.ID)
                        .fetchOne()
                        .getValue(MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCE.ID);

                for (String exonNumber : transcriptConsequence.exonNumbers()) {
                    context.insertInto(MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCEEXONNUMBER,
                            MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCEEXONNUMBER.EXONNUMBER,
                            MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCEEXONNUMBER.MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCEID)
                            .values(exonNumber, transcriptConsequenceId)
                            .execute();
                }
            }

            for (MolecularMatchParent parent : mutation.parents()) {
                int idParent = context.insertInto(MOLECULARMATCHMUTATIONPARENT,
                        MOLECULARMATCHMUTATIONPARENT.NAME,
                        MOLECULARMATCHMUTATIONPARENT.TYPE,
                        MOLECULARMATCHMUTATIONPARENT.ACTIONABLEPARENT,
                        MOLECULARMATCHMUTATIONPARENT.MOLECULARMATCHMUTATIONID)
                        .values(parent.name(), parent.type(), parent.actionableParent(), mutationId)
                        .returning(MOLECULARMATCHMUTATIONPARENT.ID)
                        .fetchOne()
                        .getValue(MOLECULARMATCHMUTATIONPARENT.ID);

                for (String transcript : parent.transcripts()) {
                    context.insertInto(MOLECULARMATCHMUTATIONPARENTTRANSCRIPT,
                            MOLECULARMATCHMUTATIONPARENTTRANSCRIPT.TRANSCRIPT,
                            MOLECULARMATCHMUTATIONPARENTTRANSCRIPT.MOLECULARMATCHMUTATIONPARENTID).values(transcript, idParent).execute();
                }
            }

            insertMutationWGSALocations(context, mutation.wgsaLocations(), mutationId);
            insertMutationWGSAMaps(context, mutation.wgsaMaps(), mutationId);
            insertMutationGRCh37Locations(context, mutation.grch37Locations(), mutationId);
            insertMutationFusions(context, mutation.fusionData(), mutationId);

            MolecularMatchExonsInfo exonsInfo = mutation.exonsInfo();
            if (exonsInfo != null) {
                insertMutationExonsInfo(context, exonsInfo, mutationId);
            }
        }
    }

    private static void insertMutationWGSALocations(@NotNull DSLContext context, @NotNull List<MolecularMatchWGSALocation> wgsaLocations,
            int mutationId) {
        for (MolecularMatchWGSALocation wgsaLocation : wgsaLocations) {
            int wgsaLocationId = context.insertInto(MOLECULARMATCHMUTATIONWGSALOCATION,
                    MOLECULARMATCHMUTATIONWGSALOCATION.CHR,
                    MOLECULARMATCHMUTATIONWGSALOCATION.START,
                    MOLECULARMATCHMUTATIONWGSALOCATION.END,
                    MOLECULARMATCHMUTATIONWGSALOCATION.REF,
                    MOLECULARMATCHMUTATIONWGSALOCATION.ALT,
                    MOLECULARMATCHMUTATIONWGSALOCATION.CHRSTARTREFALT,
                    MOLECULARMATCHMUTATIONWGSALOCATION.TRANSCRIPT,
                    MOLECULARMATCHMUTATIONWGSALOCATION.NUCLEOTIDECHANGE,
                    MOLECULARMATCHMUTATIONWGSALOCATION.AA,
                    MOLECULARMATCHMUTATIONWGSALOCATION.EXONICFUNC,
                    MOLECULARMATCHMUTATIONWGSALOCATION.POPFREQMAX,
                    MOLECULARMATCHMUTATIONWGSALOCATION.EXACAFR,
                    MOLECULARMATCHMUTATIONWGSALOCATION.EXACAMR,
                    MOLECULARMATCHMUTATIONWGSALOCATION.EXACEAS,
                    MOLECULARMATCHMUTATIONWGSALOCATION.EXACFIN,
                    MOLECULARMATCHMUTATIONWGSALOCATION.EXACNFE,
                    MOLECULARMATCHMUTATIONWGSALOCATION.EXACSAS,
                    MOLECULARMATCHMUTATIONWGSALOCATION.EXACFREQ,
                    MOLECULARMATCHMUTATIONWGSALOCATION.G1000AFR,
                    MOLECULARMATCHMUTATIONWGSALOCATION.G1000AMR,
                    MOLECULARMATCHMUTATIONWGSALOCATION.G1000EAS,
                    MOLECULARMATCHMUTATIONWGSALOCATION.G1000EUR,
                    MOLECULARMATCHMUTATIONWGSALOCATION.G1000SAS,
                    MOLECULARMATCHMUTATIONWGSALOCATION.G1000ALL,
                    MOLECULARMATCHMUTATIONWGSALOCATION.FATHMM,
                    MOLECULARMATCHMUTATIONWGSALOCATION.FATHMMPRED,
                    MOLECULARMATCHMUTATIONWGSALOCATION.ESP6500SIAA,
                    MOLECULARMATCHMUTATIONWGSALOCATION.ESP6500SIEA,
                    MOLECULARMATCHMUTATIONWGSALOCATION.DBSNP,
                    MOLECULARMATCHMUTATIONWGSALOCATION.COSMICID,
                    MOLECULARMATCHMUTATIONWGSALOCATION.PHYLOP46WAYPLACENTAL,
                    MOLECULARMATCHMUTATIONWGSALOCATION.PHYLOP100WAYVERTEBRATE,
                    MOLECULARMATCHMUTATIONWGSALOCATION.SIPHY29WAYLOGODDS,
                    MOLECULARMATCHMUTATIONWGSALOCATION.GWASSNP,
                    MOLECULARMATCHMUTATIONWGSALOCATION.GWASDIS,
                    MOLECULARMATCHMUTATIONWGSALOCATION.GWASPUBMED,
                    MOLECULARMATCHMUTATIONWGSALOCATION.GERPRS,
                    MOLECULARMATCHMUTATIONWGSALOCATION.FUNC,
                    MOLECULARMATCHMUTATIONWGSALOCATION.WGRNA,
                    MOLECULARMATCHMUTATIONWGSALOCATION.TARGETSCANS,
                    MOLECULARMATCHMUTATIONWGSALOCATION.KEYVALUE,
                    MOLECULARMATCHMUTATIONWGSALOCATION.MOLECULARMATCHMUTATIONID)
                    .values(wgsaLocation.chr(),
                            wgsaLocation.start(),
                            wgsaLocation.end(),
                            wgsaLocation.ref(),
                            wgsaLocation.alt(),
                            wgsaLocation.chrStartRefAlt(),
                            wgsaLocation.transcript(),
                            wgsaLocation.nucleotideChange(),
                            wgsaLocation.aa(),
                            wgsaLocation.exonicFunc(),
                            wgsaLocation.popFreqMax(),
                            wgsaLocation.exacAFR(),
                            wgsaLocation.exacAMR(),
                            wgsaLocation.exacEAS(),
                            wgsaLocation.exacFIN(),
                            wgsaLocation.exacNFE(),
                            wgsaLocation.exacSAS(),
                            wgsaLocation.exacFreq(),
                            wgsaLocation.g1000AFR(),
                            wgsaLocation.g1000AMR(),
                            wgsaLocation.g1000EAS(),
                            wgsaLocation.g1000EUR(),
                            wgsaLocation.g1000SAS(),
                            wgsaLocation.g1000ALL(),
                            wgsaLocation.fathmm(),
                            wgsaLocation.fathmmPred(),
                            wgsaLocation.esp6500siAA(),
                            wgsaLocation.esp6500siEA(),
                            wgsaLocation.dbSNP(),
                            wgsaLocation.cosmicId(),
                            wgsaLocation.phyloP46wayPlacental(),
                            wgsaLocation.phyloP100wayVertebrate(),
                            wgsaLocation.siPhy29wayLogOdds(),
                            wgsaLocation.gwasSNP(),
                            wgsaLocation.gwasDIS(),
                            wgsaLocation.gwasPubmed(),
                            wgsaLocation.gerpRS(),
                            wgsaLocation.func(),
                            wgsaLocation.wgRna(),
                            wgsaLocation.targetScanS(),
                            wgsaLocation.key(),
                            mutationId)
                    .returning(MOLECULARMATCHMUTATIONWGSALOCATION.ID)
                    .fetchOne()
                    .getValue(MOLECULARMATCHMUTATIONWGSALOCATION.ID);

            for (String gene : wgsaLocation.genes()) {
                context.insertInto(MOLECULARMATCHMUTATIONWGSALOCATIONGENE,
                        MOLECULARMATCHMUTATIONWGSALOCATIONGENE.GENE,
                        MOLECULARMATCHMUTATIONWGSALOCATIONGENE.MOLECULARMATCHMUTATIONWGSALOCATIONID).values(gene, wgsaLocationId).execute();
            }

            for (String fullAA : wgsaLocation.fullAAs()) {
                context.insertInto(MOLECULARMATCHMUTATIONWGSALOCATIONFULLAA,
                        MOLECULARMATCHMUTATIONWGSALOCATIONFULLAA.FULLAA,
                        MOLECULARMATCHMUTATIONWGSALOCATIONFULLAA.MOLECULARMATCHMUTATIONWGSALOCATIONID)
                        .values(fullAA, wgsaLocationId)
                        .execute();
            }

            for (String clinVarDisease : wgsaLocation.clinVarDiseases()) {
                context.insertInto(MOLECULARMATCHMUTATIONWGSALOCATIONCLINVARDISEASE,
                        MOLECULARMATCHMUTATIONWGSALOCATIONCLINVARDISEASE.CLINVARDISEASE,
                        MOLECULARMATCHMUTATIONWGSALOCATIONCLINVARDISEASE.MOLECULARMATCHMUTATIONWGSALOCATIONID)
                        .values(clinVarDisease, wgsaLocationId)
                        .execute();
            }

            for (String clinVarSig : wgsaLocation.clinVarSigs()) {
                context.insertInto(MOLECULARMATCHMUTATIONWGSALOCATIONCLINVARSIG,
                        MOLECULARMATCHMUTATIONWGSALOCATIONCLINVARSIG.CLINVARSIG,
                        MOLECULARMATCHMUTATIONWGSALOCATIONCLINVARSIG.MOLECULARMATCHMUTATIONWGSALOCATIONID)
                        .values(clinVarSig, wgsaLocationId)
                        .execute();
            }

            for (String clinVarStatus : wgsaLocation.clinVarStates()) {
                context.insertInto(MOLECULARMATCHMUTATIONWGSALOCATIONCLINVARSTATUS,
                        MOLECULARMATCHMUTATIONWGSALOCATIONCLINVARSTATUS.CLINVARSTATUS,
                        MOLECULARMATCHMUTATIONWGSALOCATIONCLINVARSTATUS.MOLECULARMATCHMUTATIONWGSALOCATIONID)
                        .values(clinVarStatus, wgsaLocationId)
                        .execute();
            }

            for (String clinVarDbId : wgsaLocation.clinVarDbIds()) {
                context.insertInto(MOLECULARMATCHMUTATIONWGSALOCATIONCLINVARDBID,
                        MOLECULARMATCHMUTATIONWGSALOCATIONCLINVARDBID.CLINVARDBID,
                        MOLECULARMATCHMUTATIONWGSALOCATIONCLINVARDBID.MOLECULARMATCHMUTATIONWGSALOCATIONID)
                        .values(clinVarDbId, wgsaLocationId)
                        .execute();
            }
        }
    }

    private static void insertMutationWGSAMaps(@NotNull DSLContext context, @NotNull List<MolecularMatchWGSAMap> wgsaMaps, int mutationId) {
        for (MolecularMatchWGSAMap wgsaMap : wgsaMaps) {
            int wgsaMapId = context.insertInto(MOLECULARMATCHMUTATIONWGSAMAP,
                    MOLECULARMATCHMUTATIONWGSAMAP.NAME,
                    MOLECULARMATCHMUTATIONWGSAMAP.GENE,
                    MOLECULARMATCHMUTATIONWGSAMAP.TRANSCRIPT,
                    MOLECULARMATCHMUTATIONWGSAMAP.EXON,
                    MOLECULARMATCHMUTATIONWGSAMAP.GRCH37CHRSTARTREFALT,
                    MOLECULARMATCHMUTATIONWGSAMAP.NUCLEOTIDECHANGE,
                    MOLECULARMATCHMUTATIONWGSAMAP.AA,
                    MOLECULARMATCHMUTATIONWGSAMAP.MOLECULARMATCHMUTATIONID)
                    .values(wgsaMap.name(),
                            wgsaMap.gene(),
                            wgsaMap.transcript(),
                            wgsaMap.exon(),
                            wgsaMap.grch37ChrStartRefAlt(),
                            wgsaMap.nucleotideChange(),
                            wgsaMap.aa(),
                            mutationId)
                    .returning(MOLECULARMATCHMUTATIONWGSAMAP.ID)
                    .fetchOne()
                    .getValue(MOLECULARMATCHMUTATIONWGSAMAP.ID);

            for (String synonym : wgsaMap.synonyms()) {
                context.insertInto(MOLECULARMATCHMUTATIONWGSAMAPSYNONYM,
                        MOLECULARMATCHMUTATIONWGSAMAPSYNONYM.SYNONYM,
                        MOLECULARMATCHMUTATIONWGSAMAPSYNONYM.MOLECULARMATCHMUTATIONWGSAMAPID).values(synonym, wgsaMapId).execute();
            }

            for (String protCoord : wgsaMap.protCoords()) {
                context.insertInto(MOLECULARMATCHMUTATIONWGSAMAPPROTCOORD,
                        MOLECULARMATCHMUTATIONWGSAMAPPROTCOORD.PROTCOORD,
                        MOLECULARMATCHMUTATIONWGSAMAPPROTCOORD.MOLECULARMATCHMUTATIONWGSAMAPID).values(protCoord, wgsaMapId).execute();
            }
        }
    }

    private static void insertMutationGRCh37Locations(@NotNull DSLContext context,
            @NotNull List<MolecularMatchGRCh37Location> grch37Locations, int mutationId) {
        for (MolecularMatchGRCh37Location grch37Location : grch37Locations) {
            int locationId = context.insertInto(MOLECULARMATCHMUTATIONGRCH37LOC,
                    MOLECULARMATCHMUTATIONGRCH37LOC.CHR,
                    MOLECULARMATCHMUTATIONGRCH37LOC.START,
                    MOLECULARMATCHMUTATIONGRCH37LOC.STOP,
                    MOLECULARMATCHMUTATIONGRCH37LOC.REF,
                    MOLECULARMATCHMUTATIONGRCH37LOC.ALT,
                    MOLECULARMATCHMUTATIONGRCH37LOC.STRAND,
                    MOLECULARMATCHMUTATIONGRCH37LOC.VALIDATED,
                    MOLECULARMATCHMUTATIONGRCH37LOC.COMPOSITEKEY,
                    MOLECULARMATCHMUTATIONGRCH37LOC.MOLECULARMATCHMUTATIONID)
                    .values(grch37Location.chr(),
                            grch37Location.start(),
                            grch37Location.stop(),
                            grch37Location.ref(),
                            grch37Location.alt(),
                            grch37Location.strand(),
                            grch37Location.validated(),
                            grch37Location.compositeKey(),
                            mutationId)
                    .returning(MOLECULARMATCHMUTATIONGRCH37LOC.ID)
                    .fetchOne()
                    .getValue(MOLECULARMATCHMUTATIONGRCH37LOC.ID);

            for (MolecularMatchGRCh37TranscriptConsequence transcriptConsequence : grch37Location.transcriptConsequences()) {
                int consequenceId = context.insertInto(MOLECULARMATCHMUTATIONGRCH37LOCCONSEQUENCE,
                        MOLECULARMATCHMUTATIONGRCH37LOCCONSEQUENCE.TRANSCRIPT,
                        MOLECULARMATCHMUTATIONGRCH37LOCCONSEQUENCE.CDNA,
                        MOLECULARMATCHMUTATIONGRCH37LOCCONSEQUENCE.AMINOACIDCHANGE,
                        MOLECULARMATCHMUTATIONGRCH37LOCCONSEQUENCE.INTRONNUMBER,
                        MOLECULARMATCHMUTATIONGRCH37LOCCONSEQUENCE.MOLECULARMATCHMUTATIONGRCH37LOCID)
                        .values(transcriptConsequence.transcript(),
                                transcriptConsequence.cdna(),
                                transcriptConsequence.aminoAcidChange(),
                                transcriptConsequence.intronNumber(),
                                locationId)
                        .returning(MOLECULARMATCHMUTATIONGRCH37LOCCONSEQUENCE.ID)
                        .fetchOne()
                        .getValue(MOLECULARMATCHMUTATIONGRCH37LOCCONSEQUENCE.ID);

                for (String txSite : transcriptConsequence.txSites()) {
                    context.insertInto(MOLECULARMATCHMUTATIONGRCH37LOCCONSEQUENCETXSITE,
                            MOLECULARMATCHMUTATIONGRCH37LOCCONSEQUENCETXSITE.TXSITE,
                            MOLECULARMATCHMUTATIONGRCH37LOCCONSEQUENCETXSITE.MOLECULARMATCHMUTATIONGRCH37LOCCONSEQUENCEID)
                            .values(txSite, consequenceId)
                            .execute();
                }

                for (String exonNumber : transcriptConsequence.exonNumbers()) {
                    context.insertInto(MOLECULARMATCHMUTATIONGRCH37LOCCONSEQUENCEEXONNUMBER,
                            MOLECULARMATCHMUTATIONGRCH37LOCCONSEQUENCEEXONNUMBER.EXONNUMBER,
                            MOLECULARMATCHMUTATIONGRCH37LOCCONSEQUENCEEXONNUMBER.MOLECULARMATCHMUTATIONGRCH37LOCCONSEQUENCEID)
                            .values(exonNumber, consequenceId)
                            .execute();
                }
            }
        }
    }

    private static void insertMutationFusions(@NotNull DSLContext context, @NotNull List<MolecularMatchFusionData> fusions,
            int mutationId) {
        for (MolecularMatchFusionData fusion : fusions) {
            int fusionId = context.insertInto(MOLECULARMATCHMUTATIONFUSION,
                    MOLECULARMATCHMUTATIONFUSION.SOURCE,
                    MOLECULARMATCHMUTATIONFUSION.SYNONYM,
                    MOLECULARMATCHMUTATIONFUSION.PAPER,
                    MOLECULARMATCHMUTATIONFUSION.MOLECULARMATCHMUTATIONID)
                    .values(fusion.source(), fusion.synonym(), fusion.paper(), mutationId)
                    .returning(MOLECULARMATCHMUTATIONFUSION.ID)
                    .fetchOne()
                    .getValue(MOLECULARMATCHMUTATIONFUSION.ID);

            for (String chromosome : fusion.aChromosomes()) {
                context.insertInto(MOLECULARMATCHMUTATIONFUSIONACHROMOSOME,
                        MOLECULARMATCHMUTATIONFUSIONACHROMOSOME.CHROMOSOME,
                        MOLECULARMATCHMUTATIONFUSIONACHROMOSOME.MOLECULARMATCHMUTATIONFUSIONID).values(chromosome, fusionId).execute();
            }

            for (String band : fusion.aBands()) {
                context.insertInto(MOLECULARMATCHMUTATIONFUSIONABAND,
                        MOLECULARMATCHMUTATIONFUSIONABAND.BAND,
                        MOLECULARMATCHMUTATIONFUSIONABAND.MOLECULARMATCHMUTATIONFUSIONID).values(band, fusionId).execute();
            }

            for (String gene : fusion.aGenes()) {
                context.insertInto(MOLECULARMATCHMUTATIONFUSIONAGENE,
                        MOLECULARMATCHMUTATIONFUSIONAGENE.GENE,
                        MOLECULARMATCHMUTATIONFUSIONAGENE.MOLECULARMATCHMUTATIONFUSIONID).values(gene, fusionId).execute();
            }

            for (String coord : fusion.aCoords()) {
                context.insertInto(MOLECULARMATCHMUTATIONFUSIONACOORD,
                        MOLECULARMATCHMUTATIONFUSIONACOORD.COORD,
                        MOLECULARMATCHMUTATIONFUSIONACOORD.MOLECULARMATCHMUTATIONFUSIONID).values(coord, fusionId).execute();
            }

            for (String transcript : fusion.aTranscripts()) {
                context.insertInto(MOLECULARMATCHMUTATIONFUSIONATRANSCRIPT,
                        MOLECULARMATCHMUTATIONFUSIONATRANSCRIPT.TRANSCRIPT,
                        MOLECULARMATCHMUTATIONFUSIONATRANSCRIPT.MOLECULARMATCHMUTATIONFUSIONID).values(transcript, fusionId).execute();
            }

            for (String orientation : fusion.aOrientations()) {
                context.insertInto(MOLECULARMATCHMUTATIONFUSIONAORIENTATION,
                        MOLECULARMATCHMUTATIONFUSIONAORIENTATION.ORIENTATION,
                        MOLECULARMATCHMUTATIONFUSIONAORIENTATION.MOLECULARMATCHMUTATIONFUSIONID).values(orientation, fusionId).execute();
            }

            for (MolecularMatchFusionGenomicRegion genomicRegion : fusion.aGenomicRegions()) {
                context.insertInto(MOLECULARMATCHMUTATIONFUSIONAGENOMICREGION,
                        MOLECULARMATCHMUTATIONFUSIONAGENOMICREGION.NUM,
                        MOLECULARMATCHMUTATIONFUSIONAGENOMICREGION.TYPE,
                        MOLECULARMATCHMUTATIONFUSIONAGENOMICREGION.MOLECULARMATCHMUTATIONFUSIONID)
                        .values(genomicRegion.num(), genomicRegion.type(), fusionId)
                        .execute();
            }

            for (String chromosome : fusion.bChromosomes()) {
                context.insertInto(MOLECULARMATCHMUTATIONFUSIONBCHROMOSOME,
                        MOLECULARMATCHMUTATIONFUSIONBCHROMOSOME.CHROMOSOME,
                        MOLECULARMATCHMUTATIONFUSIONBCHROMOSOME.MOLECULARMATCHMUTATIONFUSIONID).values(chromosome, fusionId).execute();
            }

            for (String band : fusion.bBands()) {
                context.insertInto(MOLECULARMATCHMUTATIONFUSIONBBAND,
                        MOLECULARMATCHMUTATIONFUSIONBBAND.BAND,
                        MOLECULARMATCHMUTATIONFUSIONBBAND.MOLECULARMATCHMUTATIONFUSIONID).values(band, fusionId).execute();
            }

            for (String gene : fusion.bGenes()) {
                context.insertInto(MOLECULARMATCHMUTATIONFUSIONBGENE,
                        MOLECULARMATCHMUTATIONFUSIONBGENE.GENE,
                        MOLECULARMATCHMUTATIONFUSIONBGENE.MOLECULARMATCHMUTATIONFUSIONID).values(gene, fusionId).execute();
            }

            for (String coord : fusion.bCoords()) {
                context.insertInto(MOLECULARMATCHMUTATIONFUSIONBCOORD,
                        MOLECULARMATCHMUTATIONFUSIONBCOORD.COORD,
                        MOLECULARMATCHMUTATIONFUSIONBCOORD.MOLECULARMATCHMUTATIONFUSIONID).values(coord, fusionId).execute();
            }

            for (String transcript : fusion.bTranscripts()) {
                context.insertInto(MOLECULARMATCHMUTATIONFUSIONBTRANSCRIPT,
                        MOLECULARMATCHMUTATIONFUSIONBTRANSCRIPT.TRANSCRIPT,
                        MOLECULARMATCHMUTATIONFUSIONBTRANSCRIPT.MOLECULARMATCHMUTATIONFUSIONID).values(transcript, fusionId).execute();
            }

            for (String orientation : fusion.bOrientations()) {
                context.insertInto(MOLECULARMATCHMUTATIONFUSIONBORIENTATION,
                        MOLECULARMATCHMUTATIONFUSIONBORIENTATION.ORIENTATION,
                        MOLECULARMATCHMUTATIONFUSIONBORIENTATION.MOLECULARMATCHMUTATIONFUSIONID).values(orientation, fusionId).execute();
            }

            for (MolecularMatchFusionGenomicRegion genomicRegion : fusion.bGenomicRegions()) {
                context.insertInto(MOLECULARMATCHMUTATIONFUSIONBGENOMICREGION,
                        MOLECULARMATCHMUTATIONFUSIONBGENOMICREGION.NUM,
                        MOLECULARMATCHMUTATIONFUSIONBGENOMICREGION.TYPE,
                        MOLECULARMATCHMUTATIONFUSIONBGENOMICREGION.MOLECULARMATCHMUTATIONFUSIONID)
                        .values(genomicRegion.num(), genomicRegion.type(), fusionId)
                        .execute();
            }

            for (String insert : fusion.inserts()) {
                context.insertInto(MOLECULARMATCHMUTATIONFUSIONINSERT,
                        MOLECULARMATCHMUTATIONFUSIONINSERT.INS,
                        MOLECULARMATCHMUTATIONFUSIONINSERT.MOLECULARMATCHMUTATIONFUSIONID).values(insert, fusionId).execute();
            }
        }
    }

    private static void insertMutationExonsInfo(@NotNull DSLContext context, @NotNull MolecularMatchExonsInfo exonsInfo, int mutationId) {
        int exonInfoId = context.insertInto(MOLECULARMATCHMUTATIONEXONSINFO,
                MOLECULARMATCHMUTATIONEXONSINFO.CHR,
                MOLECULARMATCHMUTATIONEXONSINFO.TRANSCRIPT,
                MOLECULARMATCHMUTATIONEXONSINFO.TXSTART,
                MOLECULARMATCHMUTATIONEXONSINFO.TXEND,
                MOLECULARMATCHMUTATIONEXONSINFO.CDSSTART,
                MOLECULARMATCHMUTATIONEXONSINFO.CDSEND,
                MOLECULARMATCHMUTATIONEXONSINFO.MOLECULARMATCHMUTATIONID)
                .values(exonsInfo.chr(),
                        exonsInfo.transcript(),
                        exonsInfo.txStart(),
                        exonsInfo.txEnd(),
                        exonsInfo.cdsStart(),
                        exonsInfo.cdsEnd(),
                        mutationId)
                .returning(MOLECULARMATCHMUTATIONEXONSINFO.ID)
                .fetchOne()
                .getValue(MOLECULARMATCHMUTATIONEXONSINFO.ID);

        MolecularMatchPosition exon1Position = exonsInfo.exonBoundaries().exon1();
        if (exon1Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON1,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON1.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON1.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON1.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon1Position.start(), exon1Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon2Position = exonsInfo.exonBoundaries().exon2();
        if (exon2Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON2,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON2.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON2.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON2.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon2Position.start(), exon2Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon3Position = exonsInfo.exonBoundaries().exon3();
        if (exon3Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON3,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON3.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON3.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON3.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon3Position.start(), exon3Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon4Position = exonsInfo.exonBoundaries().exon4();
        if (exon4Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON4,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON4.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON4.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON4.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon4Position.start(), exon4Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon5Position = exonsInfo.exonBoundaries().exon5();
        if (exon5Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON5,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON5.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON5.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON5.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon5Position.start(), exon5Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon6Position = exonsInfo.exonBoundaries().exon6();
        if (exon6Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON6,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON6.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON6.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON6.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon6Position.start(), exon6Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon7Position = exonsInfo.exonBoundaries().exon7();
        if (exon7Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON7,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON7.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON7.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON7.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon7Position.start(), exon7Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon8Position = exonsInfo.exonBoundaries().exon8();
        if (exon8Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON8,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON8.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON8.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON8.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon8Position.start(), exon8Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon9Position = exonsInfo.exonBoundaries().exon9();
        if (exon9Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON9,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON9.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON9.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON9.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon9Position.start(), exon9Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon10Position = exonsInfo.exonBoundaries().exon10();
        if (exon10Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON10,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON10.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON10.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON10.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon10Position.start(), exon10Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon11Position = exonsInfo.exonBoundaries().exon11();
        if (exon11Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON11,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON11.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON11.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON11.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon11Position.start(), exon11Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon12Position = exonsInfo.exonBoundaries().exon12();
        if (exon12Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON12,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON12.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON12.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON12.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon12Position.start(), exon12Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon13Position = exonsInfo.exonBoundaries().exon13();
        if (exon13Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON13,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON13.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON13.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON13.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon13Position.start(), exon13Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon14Position = exonsInfo.exonBoundaries().exon14();
        if (exon14Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON14,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON14.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON14.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON14.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon14Position.start(), exon14Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon15Position = exonsInfo.exonBoundaries().exon15();
        if (exon15Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON15,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON15.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON15.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON15.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon15Position.start(), exon15Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon16Position = exonsInfo.exonBoundaries().exon16();
        if (exon16Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON16,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON16.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON16.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON16.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon16Position.start(), exon16Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon17Position = exonsInfo.exonBoundaries().exon17();
        if (exon17Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON17,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON17.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON17.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON17.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon17Position.start(), exon17Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon18Position = exonsInfo.exonBoundaries().exon18();
        if (exon18Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON18,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON18.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON18.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON18.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon18Position.start(), exon18Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon19Position = exonsInfo.exonBoundaries().exon19();
        if (exon19Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON19,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON19.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON19.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON19.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon19Position.start(), exon19Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon20Position = exonsInfo.exonBoundaries().exon20();
        if (exon20Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON20,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON20.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON20.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON20.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon20Position.start(), exon20Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon21Position = exonsInfo.exonBoundaries().exon21();
        if (exon21Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON21,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON21.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON21.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON21.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon21Position.start(), exon21Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon22Position = exonsInfo.exonBoundaries().exon22();
        if (exon22Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON22,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON22.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON22.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON22.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon22Position.start(), exon22Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon23Position = exonsInfo.exonBoundaries().exon23();
        if (exon23Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON23,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON23.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON23.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON23.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon23Position.start(), exon23Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon24Position = exonsInfo.exonBoundaries().exon24();
        if (exon24Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON24,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON24.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON24.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON24.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon24Position.start(), exon24Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon25Position = exonsInfo.exonBoundaries().exon25();
        if (exon25Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON25,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON25.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON25.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON25.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon25Position.start(), exon25Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon26Position = exonsInfo.exonBoundaries().exon26();
        if (exon26Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON26,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON26.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON26.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON26.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon26Position.start(), exon26Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon27Position = exonsInfo.exonBoundaries().exon27();
        if (exon27Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON27,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON27.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON27.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON27.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon27Position.start(), exon27Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon28Position = exonsInfo.exonBoundaries().exon28();
        if (exon28Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON28,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON28.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON28.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON28.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon28Position.start(), exon28Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon29Position = exonsInfo.exonBoundaries().exon29();
        if (exon29Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON29,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON29.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON29.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON29.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon29Position.start(), exon29Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon30Position = exonsInfo.exonBoundaries().exon30();
        if (exon30Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON30,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON30.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON30.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON30.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon30Position.start(), exon30Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon31Position = exonsInfo.exonBoundaries().exon31();
        if (exon31Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON31,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON31.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON31.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON31.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon31Position.start(), exon31Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon32Position = exonsInfo.exonBoundaries().exon32();
        if (exon32Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON32,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON32.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON32.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON32.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon32Position.start(), exon32Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon33Position = exonsInfo.exonBoundaries().exon33();
        if (exon33Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON33,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON33.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON33.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON33.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon33Position.start(), exon33Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon34Position = exonsInfo.exonBoundaries().exon34();
        if (exon34Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON34,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON34.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON34.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON34.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon34Position.start(), exon34Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon35Position = exonsInfo.exonBoundaries().exon35();
        if (exon35Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON35,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON35.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON35.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON35.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon35Position.start(), exon35Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon36Position = exonsInfo.exonBoundaries().exon36();
        if (exon36Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON36,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON36.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON36.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON36.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon36Position.start(), exon36Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon37Position = exonsInfo.exonBoundaries().exon37();
        if (exon37Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON37,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON37.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON37.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON37.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon37Position.start(), exon37Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon38Position = exonsInfo.exonBoundaries().exon38();
        if (exon38Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON38,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON38.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON38.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON38.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon38Position.start(), exon38Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon39Position = exonsInfo.exonBoundaries().exon39();
        if (exon39Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON39,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON39.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON39.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON39.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon39Position.start(), exon39Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon40Position = exonsInfo.exonBoundaries().exon40();
        if (exon40Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON40,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON40.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON40.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON40.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon40Position.start(), exon40Position.stop(), exonInfoId)
                    .execute();
        }

        MolecularMatchPosition exon41Position = exonsInfo.exonBoundaries().exon41();
        if (exon41Position != null) {
            context.insertInto(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON41,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON41.START,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON41.END,
                    MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON41.MOLECULARMATCHMUTATIONEXONSINFOID)
                    .values(exon41Position.start(), exon41Position.stop(), exonInfoId)
                    .execute();
        }
    }

    static void deleteAll(@NotNull DSLContext context) {
        // Delete the exon boundaries, part of a mutation
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON1).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON2).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON3).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON4).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON5).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON6).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON7).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON8).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON9).execute();

        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON10).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON11).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON12).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON13).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON14).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON15).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON16).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON17).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON18).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON19).execute();

        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON20).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON21).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON22).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON23).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON24).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON25).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON26).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON27).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON28).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON29).execute();

        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON30).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON31).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON32).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON33).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON34).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON35).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON36).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON37).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON38).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON39).execute();

        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON40).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFOBOUNDARYEXON41).execute();

        context.deleteFrom(MOLECULARMATCHMUTATIONEXONSINFO).execute();

        // Delete the fusion info, part of mutation
        context.deleteFrom(MOLECULARMATCHMUTATIONFUSIONACHROMOSOME).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONFUSIONABAND).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONFUSIONAGENE).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONFUSIONACOORD).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONFUSIONATRANSCRIPT).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONFUSIONAORIENTATION).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONFUSIONAGENOMICREGION).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONFUSIONBCHROMOSOME).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONFUSIONBBAND).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONFUSIONBGENE).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONFUSIONBCOORD).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONFUSIONBTRANSCRIPT).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONFUSIONBORIENTATION).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONFUSIONBGENOMICREGION).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONFUSIONINSERT).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONFUSION).execute();

        // Delete the GRCh37 locations, part of mutation
        context.deleteFrom(MOLECULARMATCHMUTATIONGRCH37LOCCONSEQUENCEEXONNUMBER).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONGRCH37LOCCONSEQUENCETXSITE).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONGRCH37LOCCONSEQUENCE).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONGRCH37LOC).execute();

        // Delete the WGSA location and map tree, part of a mutation
        context.deleteFrom(MOLECULARMATCHMUTATIONWGSAMAPSYNONYM).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONWGSAMAPPROTCOORD).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONWGSAMAP).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONWGSALOCATIONCLINVARDISEASE).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONWGSALOCATIONCLINVARSIG).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONWGSALOCATIONCLINVARSTATUS).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONWGSALOCATIONCLINVARDBID).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONWGSALOCATIONFULLAA).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONWGSALOCATIONGENE).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONWGSALOCATION).execute();

        // Delete the mutation tree
        context.deleteFrom(MOLECULARMATCHMUTATIONMUTATIONTYPE).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSOURCE).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONSYNONYM).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONPATHOLOGY).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONCDNA).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCEEXONNUMBER).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONTRANSCRIPTCONSEQUENCE).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONPARENTTRANSCRIPT).execute();
        context.deleteFrom(MOLECULARMATCHMUTATIONPARENT).execute();
        context.deleteFrom(MOLECULARMATCHMUTATION).execute();

        // Delete the classification tree
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONTRANSCRIPT).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONCHROMOSOME).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONSTART).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONEND).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONREF).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONALT).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONNUCLEOTIDECHANGE).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONEXON).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONEXONICFUNC).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONPATHOLOGY).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONSOURCE).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONDBSNP).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONCOSMICID).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONPOPFREQMAX).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONPARENTTRANSCRIPT).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATIONPARENT).execute();
        context.deleteFrom(MOLECULARMATCHCLASSIFICATION).execute();

        // Delete the variant info
        context.deleteFrom(MOLECULARMATCHVARIANTINFOLOCATIONEXONNUMBER).execute();
        context.deleteFrom(MOLECULARMATCHVARIANTINFOLOCATION).execute();
        context.deleteFrom(MOLECULARMATCHVARIANTINFOCONSEQUENCE).execute();
        context.deleteFrom(MOLECULARMATCHVARIANTINFOFUSION).execute();
        context.deleteFrom(MOLECULARMATCHVARIANTINFO).execute();

        // Delete the complex tables directly dependent on molecular match
        context.deleteFrom(MOLECULARMATCHSOURCE).execute();
        context.deleteFrom(MOLECULARMATCHTIEREXPLANATION).execute();
        context.deleteFrom(MOLECULARMATCHTHERAPEUTICCONTEXT).execute();
        context.deleteFrom(MOLECULARMATCHTAG).execute();
        context.deleteFrom(MOLECULARMATCHCRITERIAUNMET).execute();
        context.deleteFrom(MOLECULARMATCHPREVALENCE).execute();

        // Delete the flat tables directly dependent on molecular match
        context.deleteFrom(MOLECULARMATCHINCLUDEGENE1).execute();
        context.deleteFrom(MOLECULARMATCHINCLUDEFINDING1).execute();
        context.deleteFrom(MOLECULARMATCHINCLUDECONDITION1).execute();
        context.deleteFrom(MOLECULARMATCHINCLUDEMUTATION1).execute();
        context.deleteFrom(MOLECULARMATCHINCLUDEDRUG1).execute();
        context.deleteFrom(MOLECULARMATCHINCLUDEDRUGCLASS1).execute();
        context.deleteFrom(MOLECULARMATCHINCLUDERESISTANCE1).execute();
        context.deleteFrom(MOLECULARMATCHINCLUDESTAGE0).execute();
        context.deleteFrom(MOLECULARMATCHINCLUDEGENE0).execute();
        context.deleteFrom(MOLECULARMATCHINCLUDECONDITION0).execute();
        context.deleteFrom(MOLECULARMATCHINCLUDEMUTATION0).execute();
        context.deleteFrom(MOLECULARMATCHCRITERIAMET).execute();
        context.deleteFrom(MOLECULARMATCHEXTERNALID).execute();

        // Delete the ast tree
        context.deleteFrom(MOLECULARMATCHASTRIGHTLEFT).execute();
        context.deleteFrom(MOLECULARMATCHASTRIGHTRIGHT).execute();
        context.deleteFrom(MOLECULARMATCHASTRIGHT).execute();
        context.deleteFrom(MOLECULARMATCHASTLEFTLEFT).execute();
        context.deleteFrom(MOLECULARMATCHASTLEFTRIGHT).execute();
        context.deleteFrom(MOLECULARMATCHASTLEFT).execute();
        context.deleteFrom(MOLECULARMATCHAST).execute();

        // Finally delete the main node
        context.deleteFrom(MOLECULARMATCH).execute();
    }
}
