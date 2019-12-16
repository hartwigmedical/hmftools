package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.BRCA;
import static com.hartwig.hmftools.vicc.database.Tables.BRCAANNOTATION1000GENOMES;
import static com.hartwig.hmftools.vicc.database.Tables.BRCAANNOTATIONBIC;
import static com.hartwig.hmftools.vicc.database.Tables.BRCAANNOTATIONCLINVAR;
import static com.hartwig.hmftools.vicc.database.Tables.BRCAANNOTATIONENIGMA;
import static com.hartwig.hmftools.vicc.database.Tables.BRCAANNOTATIONESP;
import static com.hartwig.hmftools.vicc.database.Tables.BRCAANNOTATIONEXAC;
import static com.hartwig.hmftools.vicc.database.Tables.BRCAANNOTATIONEXLOVD;
import static com.hartwig.hmftools.vicc.database.Tables.BRCAANNOTATIONLOVD;

import com.hartwig.hmftools.vicc.datamodel.brca.Brca;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

final class BRCADAOFunctions {

    private BRCADAOFunctions() {
    }

    static void write(@NotNull DSLContext context, int viccEntryId, @NotNull Brca brca) {
        int id = context.insertInto(BRCA,
                BRCA.GENESYMBOL,
                BRCA.CHR,
                BRCA.POS,
                BRCA.REF,
                BRCA.ALT,
                BRCA.GENOMICCOORDINATEHG36,
                BRCA.HG36START,
                BRCA.HG36END,
                BRCA.GENOMICCOORDINATEHG37,
                BRCA.HG37START,
                BRCA.HG37END,
                BRCA.GENOMICCOORDINATEHG38,
                BRCA.HG38START,
                BRCA.HG38END,
                BRCA.PROTEINCHANGE,
                BRCA.REFERENCESEQUENCE,
                BRCA.SYNONYMS,
                BRCA.HGVSCDNA,
                BRCA.HGVSPROTEIN,
                BRCA.HGVSRNA,
                BRCA.SIFTSCORE,
                BRCA.SIFTPREDICTION,
                BRCA.POLYPHENSCORE,
                BRCA.POLYPHENPREDICTION,
                BRCA.PATHOGENICITYALL,
                BRCA.PATHOGENICITYEXPERT,
                BRCA.ALLELEFREQUENCY,
                BRCA.MAXALLELEFREQUENCY,
                BRCA.DISCORDANT,
                BRCA.IDBRCA,
                BRCA.CHANGETYPEID,
                BRCA.DATARELEASEID,
                BRCA.SOURCE,
                BRCA.SOURCEURL,
                BRCA.VICCENTRYID)
                .values(brca.geneSymbol(),
                        brca.chr(),
                        brca.pos(),
                        brca.ref(),
                        brca.alt(),
                        brca.genomicCoordinateHg36(),
                        brca.hg36Start(),
                        brca.hg36End(),
                        brca.genomicCoordinateHg37(),
                        brca.hg37Start(),
                        brca.hg37End(),
                        brca.genomicCoordinateHg38(),
                        brca.hg38Start(),
                        brca.hg38End(),
                        brca.proteinChange(),
                        brca.referenceSequence(),
                        brca.synonyms(),
                        brca.hgvsCDNA(),
                        brca.hgvsProtein(),
                        brca.hgvsRNA(),
                        brca.siftScore(),
                        brca.siftPrediction(),
                        brca.polyphenScore(),
                        brca.polyphenPrediction(),
                        brca.pathogenicityAll(),
                        brca.pathogenicityExpert(),
                        brca.alleleFrequency(),
                        brca.maxAlleleFrequency(),
                        brca.discordant(),
                        brca.id(),
                        brca.changeTypeId(),
                        brca.dataReleaseId(),
                        brca.source(),
                        brca.sourceURL(),
                        viccEntryId)
                .returning(BRCA.ID)
                .fetchOne()
                .getValue(BRCA.ID);

        context.insertInto(BRCAANNOTATION1000GENOMES,
                BRCAANNOTATION1000GENOMES.VARIANTIN1000GENOMES,
                BRCAANNOTATION1000GENOMES.BXID,
                BRCAANNOTATION1000GENOMES.ALLELEFREQUENCY,
                BRCAANNOTATION1000GENOMES.AFRALLELEFREQUENCY,
                BRCAANNOTATION1000GENOMES.AMRALLELEFREQUENCY,
                BRCAANNOTATION1000GENOMES.EASALLELEFREQUENCY,
                BRCAANNOTATION1000GENOMES.EURALLELEFREQUENCY,
                BRCAANNOTATION1000GENOMES.SASALLELEFREQUENCY,
                BRCAANNOTATION1000GENOMES.BRCAID)
                .values(brca.annotation1000Genomes().variantIn1000Genomes(),
                        brca.annotation1000Genomes().bxId(),
                        brca.annotation1000Genomes().alleleFrequency(),
                        brca.annotation1000Genomes().afrAlleleFrequency(),
                        brca.annotation1000Genomes().amrAlleleFrequency(),
                        brca.annotation1000Genomes().easAlleleFrequency(),
                        brca.annotation1000Genomes().eurAlleleFrequency(),
                        brca.annotation1000Genomes().sasAlleleFrequency(),
                        id)
                .execute();

        context.insertInto(BRCAANNOTATIONBIC,
                BRCAANNOTATIONBIC.VARIANTINBIC,
                BRCAANNOTATIONBIC.BXID,
                BRCAANNOTATIONBIC.MUTATIONTYPE,
                BRCAANNOTATIONBIC.CLINICALCLASSIFICATION,
                BRCAANNOTATIONBIC.CLINICALIMPORTANCE,
                BRCAANNOTATIONBIC.NOMENCLATURE,
                BRCAANNOTATIONBIC.ETHNICITY,
                BRCAANNOTATIONBIC.PATIENTNATIONALITY,
                BRCAANNOTATIONBIC.GERMLINEORSOMATIC,
                BRCAANNOTATIONBIC.NUMBEROFFAMILYMEMBERCARRYINGMUTATION,
                BRCAANNOTATIONBIC.LITERATURECITATION,
                BRCAANNOTATIONBIC.BRCAID)
                .values(brca.annotationBIC().variantInBIC(),
                        brca.annotationBIC().bxId(),
                        brca.annotationBIC().mutationType(),
                        brca.annotationBIC().clinicalClassification(),
                        brca.annotationBIC().clinicalImportance(),
                        brca.annotationBIC().nomenclature(),
                        brca.annotationBIC().ethnicity(),
                        brca.annotationBIC().patientNationality(),
                        brca.annotationBIC().germlineOrSomatic(),
                        brca.annotationBIC().numberOfFamilyMemberCarryingMutation(),
                        brca.annotationBIC().literatureCitation(),
                        id)
                .execute();

        context.insertInto(BRCAANNOTATIONCLINVAR,
                BRCAANNOTATIONCLINVAR.VARIANTINCLINVAR,
                BRCAANNOTATIONCLINVAR.BXID,
                BRCAANNOTATIONCLINVAR.CLINICALSIGNIFICANCE,
                BRCAANNOTATIONCLINVAR.SUBMITTER,
                BRCAANNOTATIONCLINVAR.METHOD,
                BRCAANNOTATIONCLINVAR.ALLELEORIGIN,
                BRCAANNOTATIONCLINVAR.SCV,
                BRCAANNOTATIONCLINVAR.DATELASTUPDATED,
                BRCAANNOTATIONCLINVAR.BRCAID)
                .values(brca.annotationClinVar().variantInClinVar(),
                        brca.annotationClinVar().bxId(),
                        brca.annotationClinVar().clinicalSignificance(),
                        brca.annotationClinVar().submitter(),
                        brca.annotationClinVar().method(),
                        brca.annotationClinVar().alleleOrigin(),
                        brca.annotationClinVar().scv(),
                        brca.annotationClinVar().dateLastUpdated(),
                        id)
                .execute();

        context.insertInto(BRCAANNOTATIONENIGMA,
                BRCAANNOTATIONENIGMA.VARIANTINENIGMA,
                BRCAANNOTATIONENIGMA.BXID,
                BRCAANNOTATIONENIGMA.ALLELEORIGIN,
                BRCAANNOTATIONENIGMA.CLINVARACCESSION,
                BRCAANNOTATIONENIGMA.ASSERTIONMETHOD,
                BRCAANNOTATIONENIGMA.ASSERTIONMETHODCITATION,
                BRCAANNOTATIONENIGMA.COLLECTIONMETHOD,
                BRCAANNOTATIONENIGMA.CONDITIONCATEGORY,
                BRCAANNOTATIONENIGMA.CONDITIONIDVALUE,
                BRCAANNOTATIONENIGMA.CONDITIONIDTYPE,
                BRCAANNOTATIONENIGMA.CLINICALSIGNIFICANCE,
                BRCAANNOTATIONENIGMA.CLINICALSIGNIFICANCECITATIONS,
                BRCAANNOTATIONENIGMA.COMMENTONCLINICALSIGNIFICANCE,
                BRCAANNOTATIONENIGMA.DATELASTEVALUATED,
                BRCAANNOTATIONENIGMA.URL,
                BRCAANNOTATIONENIGMA.BRCAID)
                .values(brca.annotationENIGMA().variantInENIGMA(),
                        brca.annotationENIGMA().bxId(),
                        brca.annotationENIGMA().alleleOrigin(),
                        brca.annotationENIGMA().clinVarAccession(),
                        brca.annotationENIGMA().assertionMethod(),
                        brca.annotationENIGMA().assertionMethodCitation(),
                        brca.annotationENIGMA().collectionMethod(),
                        brca.annotationENIGMA().conditionCategory(),
                        brca.annotationENIGMA().conditionIdValue(),
                        brca.annotationENIGMA().conditionIdType(),
                        brca.annotationENIGMA().clinicalSignificance(),
                        brca.annotationENIGMA().clinicalSignificanceCitations(),
                        brca.annotationENIGMA().commentOnClinicalSignificance(),
                        brca.annotationENIGMA().dateLastEvaluated(),
                        brca.annotationENIGMA().url(),
                        id)
                .execute();

        context.insertInto(BRCAANNOTATIONESP,
                BRCAANNOTATIONESP.VARIANTINESP,
                BRCAANNOTATIONESP.BXID,
                BRCAANNOTATIONESP.MINORALLELEFREQUENCYPERCENT,
                BRCAANNOTATIONESP.ALLELEFREQUENCY,
                BRCAANNOTATIONESP.AAALLELEFREQUENCY,
                BRCAANNOTATIONESP.EAALLELEFREQUENCY,
                BRCAANNOTATIONESP.BRCAID)
                .values(brca.annotationESP().variantInESP(),
                        brca.annotationESP().bxId(),
                        brca.annotationESP().minorAlleleFrequencyPercent(),
                        brca.annotationESP().alleleFrequency(),
                        brca.annotationESP().aaAlleleFrequency(),
                        brca.annotationESP().eaAlleleFrequency(),
                        id)
                .execute();

        context.insertInto(BRCAANNOTATIONEXAC,
                BRCAANNOTATIONEXAC.VARIANTINEXAC,
                BRCAANNOTATIONEXAC.BXID,
                BRCAANNOTATIONEXAC.ALLELEFREQUENCY,
                BRCAANNOTATIONEXAC.ALLELEFREQUENCYAFR,
                BRCAANNOTATIONEXAC.ALLELEFREQUENCYAMR,
                BRCAANNOTATIONEXAC.ALLELEFREQUENCYEAS,
                BRCAANNOTATIONEXAC.ALLELEFREQUENCYFIN,
                BRCAANNOTATIONEXAC.ALLELEFREQUENCYNFE,
                BRCAANNOTATIONEXAC.ALLELEFREQUENCYOTH,
                BRCAANNOTATIONEXAC.ALLELEFREQUENCYSAS,
                BRCAANNOTATIONEXAC.ALLELENUMBERAFR,
                BRCAANNOTATIONEXAC.ALLELENUMBERAMR,
                BRCAANNOTATIONEXAC.ALLELENUMBEREAS,
                BRCAANNOTATIONEXAC.ALLELENUMBERFIN,
                BRCAANNOTATIONEXAC.ALLELENUMBERNFE,
                BRCAANNOTATIONEXAC.ALLELENUMBEROTH,
                BRCAANNOTATIONEXAC.ALLELENUMBERSAS,
                BRCAANNOTATIONEXAC.HOMOZYGOUSCOUNTAFR,
                BRCAANNOTATIONEXAC.HOMOZYGOUSCOUNTAMR,
                BRCAANNOTATIONEXAC.HOMOZYGOUSCOUNTEAS,
                BRCAANNOTATIONEXAC.HOMOZYGOUSCOUNTFIN,
                BRCAANNOTATIONEXAC.HOMOZYGOUSCOUNTNFE,
                BRCAANNOTATIONEXAC.HOMOZYGOUSCOUNTOTH,
                BRCAANNOTATIONEXAC.HOMOZYGOUSCOUNTSAS,
                BRCAANNOTATIONEXAC.ALLELECOUNTAFR,
                BRCAANNOTATIONEXAC.ALLELECOUNTAMR,
                BRCAANNOTATIONEXAC.ALLELECOUNTEAS,
                BRCAANNOTATIONEXAC.ALLELECOUNTFIN,
                BRCAANNOTATIONEXAC.ALLELECOUNTNFE,
                BRCAANNOTATIONEXAC.ALLELECOUNTOTH,
                BRCAANNOTATIONEXAC.ALLELECOUNTSAS,
                BRCAANNOTATIONEXAC.BRCAID)
                .values(brca.annotationExAC().variantInExAC(),
                        brca.annotationExAC().bxId(),
                        brca.annotationExAC().alleleFrequency(),
                        brca.annotationExAC().alleleFrequencyAFR(),
                        brca.annotationExAC().alleleFrequencyAMR(),
                        brca.annotationExAC().alleleFrequencyEAS(),
                        brca.annotationExAC().alleleFrequencyFIN(),
                        brca.annotationExAC().alleleFrequencyNFE(),
                        brca.annotationExAC().alleleFrequencyOTH(),
                        brca.annotationExAC().alleleFrequencySAS(),
                        brca.annotationExAC().alleleNumberAFR(),
                        brca.annotationExAC().alleleNumberAMR(),
                        brca.annotationExAC().alleleNumberEAS(),
                        brca.annotationExAC().alleleNumberFIN(),
                        brca.annotationExAC().alleleNumberNFE(),
                        brca.annotationExAC().alleleNumberOTH(),
                        brca.annotationExAC().alleleNumberSAS(),
                        brca.annotationExAC().homozygousCountAFR(),
                        brca.annotationExAC().homozygousCountAMR(),
                        brca.annotationExAC().homozygousCountEAS(),
                        brca.annotationExAC().homozygousCountFIN(),
                        brca.annotationExAC().homozygousCountNFE(),
                        brca.annotationExAC().homozygousCountOTH(),
                        brca.annotationExAC().homozygousCountSAS(),
                        brca.annotationExAC().alleleCountAFR(),
                        brca.annotationExAC().alleleCountAMR(),
                        brca.annotationExAC().alleleCountEAS(),
                        brca.annotationExAC().alleleCountFIN(),
                        brca.annotationExAC().alleleCountNFE(),
                        brca.annotationExAC().alleleCountOTH(),
                        brca.annotationExAC().alleleCountSAS(),
                        id)
                .execute();

        context.insertInto(BRCAANNOTATIONEXLOVD,
                BRCAANNOTATIONEXLOVD.VARIANTINEXLOVD,
                BRCAANNOTATIONEXLOVD.BXID,
                BRCAANNOTATIONEXLOVD.COOCCURRENCELR,
                BRCAANNOTATIONEXLOVD.SUMFAMILYLR,
                BRCAANNOTATIONEXLOVD.SEGREGATIONLR,
                BRCAANNOTATIONEXLOVD.POSTERIORPROBABILITY,
                BRCAANNOTATIONEXLOVD.MISSENSEANALYSISPRIORPROBABILITY,
                BRCAANNOTATIONEXLOVD.COMBINEDPRIORPROBABILITY,
                BRCAANNOTATIONEXLOVD.IARCCLASS,
                BRCAANNOTATIONEXLOVD.LITERATURESOURCE,
                BRCAANNOTATIONEXLOVD.BRCAID)
                .values(brca.annotationExLOVD().variantInExLOVD(),
                        brca.annotationExLOVD().bxId(),
                        brca.annotationExLOVD().cooccurrenceLR(),
                        brca.annotationExLOVD().sumFamilyLR(),
                        brca.annotationExLOVD().segregationLR(),
                        brca.annotationExLOVD().posteriorProbability(),
                        brca.annotationExLOVD().missenseAnalysisPriorProbability(),
                        brca.annotationExLOVD().combinedPriorProbability(),
                        brca.annotationExLOVD().iarcClass(),
                        brca.annotationExLOVD().literatureSource(),
                        id)
                .execute();

        context.insertInto(BRCAANNOTATIONLOVD,
                BRCAANNOTATIONLOVD.VARIANTINLOVD,
                BRCAANNOTATIONLOVD.BXID,
                BRCAANNOTATIONLOVD.DBID,
                BRCAANNOTATIONLOVD.HGVSCDNA,
                BRCAANNOTATIONLOVD.HGVSPROTEIN,
                BRCAANNOTATIONLOVD.RNA,
                BRCAANNOTATIONLOVD.VARIANTEFFECT,
                BRCAANNOTATIONLOVD.VARIANTFREQUENCY,
                BRCAANNOTATIONLOVD.VARIANTHAPLOTYPE,
                BRCAANNOTATIONLOVD.GENETICORIGIN,
                BRCAANNOTATIONLOVD.FUNCTIONALANALYSISTECHNIQUE,
                BRCAANNOTATIONLOVD.FUNCTIONALANALYSISRESULT,
                BRCAANNOTATIONLOVD.SUBMITTERS,
                BRCAANNOTATIONLOVD.INDIVIDUALS,
                BRCAANNOTATIONLOVD.BRCAID)
                .values(brca.annotationLOVD().variantInLOVD(),
                        brca.annotationLOVD().bxId(),
                        brca.annotationLOVD().dbId(),
                        brca.annotationLOVD().hgvsCDNA(),
                        brca.annotationLOVD().hgvsProtein(),
                        brca.annotationLOVD().rna(),
                        brca.annotationLOVD().variantEffect(),
                        brca.annotationLOVD().variantFrequency(),
                        brca.annotationLOVD().variantHaplotype(),
                        brca.annotationLOVD().geneticOrigin(),
                        brca.annotationLOVD().functionalAnalysisTechnique(),
                        brca.annotationLOVD().functionalAnalysisResult(),
                        brca.annotationLOVD().submitters(),
                        brca.annotationLOVD().individuals(),
                        id)
                .execute();
    }

    static void deleteAll(@NotNull DSLContext context) {
        // First delete the tables dependent on BRCA
        context.deleteFrom(BRCAANNOTATION1000GENOMES).execute();
        context.deleteFrom(BRCAANNOTATIONBIC).execute();
        context.deleteFrom(BRCAANNOTATIONCLINVAR).execute();
        context.deleteFrom(BRCAANNOTATIONENIGMA).execute();
        context.deleteFrom(BRCAANNOTATIONESP).execute();
        context.deleteFrom(BRCAANNOTATIONEXAC).execute();
        context.deleteFrom(BRCAANNOTATIONEXLOVD).execute();
        context.deleteFrom(BRCAANNOTATIONLOVD).execute();

        // Then delete the main object
        context.deleteFrom(BRCA).execute();
    }
}
