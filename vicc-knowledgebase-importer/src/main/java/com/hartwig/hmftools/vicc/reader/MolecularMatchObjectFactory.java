package com.hartwig.hmftools.vicc.reader;

import static com.hartwig.hmftools.common.utils.json.JsonFunctions.nullableString;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.optionalJsonArray;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.optionalJsonObject;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.optionalNullableString;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.optionalString;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.optionalStringList;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.string;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.stringList;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatch;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAst;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAstLeft;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAstLeftLeft;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAstLeftRight;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAstRight;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAstRightLeft;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAstRightRight;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchClassification;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchCriteriaUnmet;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchExonBoundaries;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchExonsInfo;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchFusion;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchFusionData;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchFusionGenomicRegion;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchGRCh37Location;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchGRCh37TranscriptConsequence;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchLocation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchMutation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchParent;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchPosition;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchPrevalence;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchSource;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchTag;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchTherapeuticContext;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchTierExplanation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchTranscriptConsequence;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchVariantInfo;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchWGSALocation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchWGSAMap;
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
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchExonBoundaries;
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
import org.jetbrains.annotations.Nullable;

final class MolecularMatchObjectFactory {

    private MolecularMatchObjectFactory() {
    }

    @NotNull
    static MolecularMatch create(@NotNull JsonObject molecularMatchObject) {
        ViccDatamodelCheckerFactory.molecularMatchEntryChecker().check(molecularMatchObject);

        return ImmutableMolecularMatch.builder()
                .direction(string(molecularMatchObject, "direction"))
                .biomarkerClass(string(molecularMatchObject, "biomarkerClass"))
                .mutations(createMutations(molecularMatchObject.getAsJsonArray("mutations")))
                .variantInfos(createVariantInfos(molecularMatchObject.getAsJsonArray("variantInfo")))
                .prevalences(createPrevalences(molecularMatchObject.getAsJsonArray("prevalence")))
                .score(string(molecularMatchObject, "_score"))
                .sources(createSources(molecularMatchObject.getAsJsonArray("sources")))
                .clinicalSignificance(string(molecularMatchObject, "clinicalSignificance"))
                .tier(string(molecularMatchObject, "tier"))
                .tierExplanations(createTierExplanations(molecularMatchObject.getAsJsonArray("tierExplanation")))
                .ampcap(string(molecularMatchObject, "ampcap"))
                .civicValue(string(molecularMatchObject, "civic"))
                .regulatoryBody(string(molecularMatchObject, "regulatoryBody"))
                .regulatoryBodyApproved(string(molecularMatchObject, "regulatoryBodyApproved"))
                .guidelineBody(optionalString(molecularMatchObject, "guidelineBody"))
                .guidelineVersion(optionalString(molecularMatchObject, "guidelineVersion"))
                .includeGene1(optionalStringList(molecularMatchObject, "includeGene1"))
                .includeFinding1(optionalStringList(molecularMatchObject, "includeFinding1"))
                .includeCondition1(stringList(molecularMatchObject, "includeCondition1"))
                .includeMutation1(optionalStringList(molecularMatchObject, "includeMutation1"))
                .includeDrug1(optionalStringList(molecularMatchObject, "includeDrug1"))
                .includeDrugClass1(optionalStringList(molecularMatchObject, "includeDrugclass1"))
                .includeResistance1(optionalStringList(molecularMatchObject, "includeResistance1"))
                .includeStage0(optionalStringList(molecularMatchObject, "includeStage0"))
                .includeGene0(optionalStringList(molecularMatchObject, "includeGene0"))
                .includeCondition0(stringList(molecularMatchObject, "includeCondition0"))
                .includeMutation0(optionalStringList(molecularMatchObject, "includeMutation0"))
                .criteriaMets(stringList(molecularMatchObject, "criteriaMet"))
                .criteriaUnmets(createCriteriaUnmets(molecularMatchObject.getAsJsonArray("criteriaUnmet")))
                .ast(createAst(molecularMatchObject.getAsJsonObject("ast")))
                .institutions(optionalStringList(molecularMatchObject, "institution"))
                .tags(createTags(molecularMatchObject.getAsJsonArray("tags")))
                .classifications(createClassifications(molecularMatchObject.getAsJsonArray("classifications")))
                .noTherapyAvailable(optionalString(molecularMatchObject, "noTherapyAvailable"))
                .therapeuticContexts(createTherapeuticContexts(molecularMatchObject.getAsJsonArray("therapeuticContext")))
                .sixtier(string(molecularMatchObject, "sixtier"))
                .mvld(string(molecularMatchObject, "mvld"))
                .autoGenerateNarrative(string(molecularMatchObject, "autoGenerateNarrative"))
                .narrative(string(molecularMatchObject, "narrative"))
                .expression(string(molecularMatchObject, "expression"))
                .customer(string(molecularMatchObject, "customer"))
                .version(string(molecularMatchObject, "version"))
                .id(string(molecularMatchObject, "id"))
                .externalIds(optionalStringList(molecularMatchObject, "external_id"))
                .uniqueKey(string(molecularMatchObject, "uniqueKey"))
                .hashKey(string(molecularMatchObject, "hashKey"))
                .build();
    }

    @NotNull
    private static List<MolecularMatchMutation> createMutations(@NotNull JsonArray mutationArray) {
        List<MolecularMatchMutation> mutationList = Lists.newArrayList();
        ViccDatamodelChecker mutationChecker = ViccDatamodelCheckerFactory.molecularMatchMutationChecker();

        for (JsonElement mutationElement : mutationArray) {
            JsonObject mutationObject = mutationElement.getAsJsonObject();
            mutationChecker.check(mutationObject);

            mutationList.add(ImmutableMolecularMatchMutation.builder()
                    .geneSymbol(string(mutationObject, "geneSymbol"))
                    .name(string(mutationObject, "name"))
                    .transcriptRecognized(optionalString(mutationObject, "transcriptRecognized"))
                    .transcript(optionalString(mutationObject, "transcript"))
                    .longestTranscript(optionalString(mutationObject, "longestTranscript"))
                    .uniprotTranscript(optionalString(mutationObject, "uniprotTranscript"))
                    .transcriptConsequences(createTranscriptConsequences(optionalJsonArray(mutationObject, "transcriptConsequence")))
                    .parents(createParents(mutationObject.getAsJsonArray("parents")))
                    .wgsaLocations(createWGSALocations(optionalJsonObject(mutationObject, "wgsaData")))
                    .wgsaMaps(createWGSAMaps(optionalJsonArray(mutationObject, "wgsaMap")))
                    .exonsInfo(createExonsInfo(optionalJsonObject(mutationObject, "exonsInfo")))
                    .fusionData(createFusionData(optionalJsonArray(mutationObject, "fusionData")))
                    .mutationTypes(stringList(mutationObject, "mutation_type"))
                    .sources(stringList(mutationObject, "sources"))
                    .synonyms(stringList(mutationObject, "synonyms"))
                    .grch37Locations(createGRCh37Locations(mutationObject.getAsJsonArray("GRCh37_location")))
                    .pathology(stringList(mutationObject, "pathology"))
                    .cDNA(stringList(mutationObject, "cdna"))
                    .description(string(mutationObject, "description"))
                    .src(string(mutationObject, "_src"))
                    .id(string(mutationObject, "id"))
                    .build());
        }
        return mutationList;
    }

    @NotNull
    private static List<MolecularMatchTranscriptConsequence> createTranscriptConsequences(@Nullable JsonArray transcriptConsequenceArray) {
        if (transcriptConsequenceArray == null) {
            return Lists.newArrayList();
        }

        List<MolecularMatchTranscriptConsequence> transcriptConsequenceList = Lists.newArrayList();
        ViccDatamodelChecker transcriptConsequenceChecker = ViccDatamodelCheckerFactory.molecularMatchTranscriptConsequenceChecker();

        for (JsonElement transcriptConsequenceElement : transcriptConsequenceArray) {
            JsonObject transcriptConsequenceObject = transcriptConsequenceElement.getAsJsonObject();
            transcriptConsequenceChecker.check(transcriptConsequenceObject);

            transcriptConsequenceList.add(ImmutableMolecularMatchTranscriptConsequence.builder()
                    .chr(optionalString(transcriptConsequenceObject, "chr"))
                    .start(optionalString(transcriptConsequenceObject, "start"))
                    .stop(optionalString(transcriptConsequenceObject, "stop"))
                    .ref(optionalString(transcriptConsequenceObject, "ref"))
                    .alt(optionalString(transcriptConsequenceObject, "alt"))
                    .referenceGenome(string(transcriptConsequenceObject, "referenceGenome"))
                    .transcript(string(transcriptConsequenceObject, "transcript"))
                    .strand(string(transcriptConsequenceObject, "strand"))
                    .cdna(optionalString(transcriptConsequenceObject, "cdna"))
                    .aminoAcidChange(optionalNullableString(transcriptConsequenceObject, "amino_acid_change"))
                    .intronNumber(nullableString(transcriptConsequenceObject, "intronNumber"))
                    .exonNumbers(optionalStringList(transcriptConsequenceObject, "exonNumber"))
                    .suppress(string(transcriptConsequenceObject, "suppress"))
                    .custom(string(transcriptConsequenceObject, "custom"))
                    .validated(string(transcriptConsequenceObject, "validated"))
                    .compositeKey(string(transcriptConsequenceObject, "compositeKey"))
                    .build());
        }
        return transcriptConsequenceList;
    }

    @NotNull
    private static List<MolecularMatchWGSALocation> createWGSALocations(@Nullable JsonObject wgsaDataObject) {
        if (wgsaDataObject == null) {
            return Lists.newArrayList();
        }

        ViccDatamodelCheckerFactory.molecularMatchWGSADataChecker().check(wgsaDataObject);

        List<MolecularMatchWGSALocation> wgsaLocationList = Lists.newArrayList();
        ViccDatamodelChecker wgsaLocationChecker = ViccDatamodelCheckerFactory.molecularMatchWGSALocationChecker();

        for (JsonElement wgsaLocationElement : wgsaDataObject.get("locations").getAsJsonArray()) {
            JsonObject wgsaLocationObject = wgsaLocationElement.getAsJsonObject();
            wgsaLocationChecker.check(wgsaLocationObject);

            wgsaLocationList.add(ImmutableMolecularMatchWGSALocation.builder()
                    .genes(stringList(wgsaLocationObject, "Gene"))
                    .chr(string(wgsaLocationObject, "Chr"))
                    .start(string(wgsaLocationObject, "Start"))
                    .end(string(wgsaLocationObject, "End"))
                    .ref(string(wgsaLocationObject, "Ref"))
                    .alt(string(wgsaLocationObject, "Alt"))
                    .chrStartRefAlt(string(wgsaLocationObject, "Chr_Start_Ref_Alt"))
                    .transcript(string(wgsaLocationObject, "Transcript"))
                    .nucleotideChange(string(wgsaLocationObject, "NucleotideChange"))
                    .aa(optionalString(wgsaLocationObject, "AA"))
                    .fullAAs(stringList(wgsaLocationObject, "FullAA"))
                    .exonicFunc(optionalString(wgsaLocationObject, "ExonicFunc"))
                    .popFreqMax(string(wgsaLocationObject, "PopFreqMax"))
                    .clinVarDiseases(optionalStringList(wgsaLocationObject, "ClinVar_DIS"))
                    .clinVarSigs(optionalStringList(wgsaLocationObject, "ClinVar_SIG"))
                    .clinVarStates(optionalStringList(wgsaLocationObject, "ClinVar_STATUS"))
                    .clinVarDbIds(optionalStringList(wgsaLocationObject, "ClinVar_DBID"))
                    .exacAFR(optionalString(wgsaLocationObject, "ExAC_AFR"))
                    .exacAMR(optionalString(wgsaLocationObject, "ExAC_AMR"))
                    .exacEAS(optionalString(wgsaLocationObject, "ExAC_EAS"))
                    .exacFIN(optionalString(wgsaLocationObject, "ExAC_FIN"))
                    .exacNFE(optionalString(wgsaLocationObject, "ExAC_NFE"))
                    .exacSAS(optionalString(wgsaLocationObject, "ExAC_SAS"))
                    .exacFreq(optionalString(wgsaLocationObject, "ExAC_Freq"))
                    .g1000AFR(optionalString(wgsaLocationObject, "1000G_AFR"))
                    .g1000AMR(optionalString(wgsaLocationObject, "1000G_AMR"))
                    .g1000EUR(optionalString(wgsaLocationObject, "1000G_EUR"))
                    .g1000EAS(optionalString(wgsaLocationObject, "1000G_EAS"))
                    .g1000SAS(optionalString(wgsaLocationObject, "1000G_SAS"))
                    .g1000ALL(optionalString(wgsaLocationObject, "1000G_ALL"))
                    .fathmm(string(wgsaLocationObject, "FATHMM"))
                    .fathmmPred(string(wgsaLocationObject, "FATHMM_Pred"))
                    .esp6500siAA(optionalString(wgsaLocationObject, "ESP6500si_AA"))
                    .esp6500siEA(optionalString(wgsaLocationObject, "ESP6500si_EA"))
                    .dbSNP(optionalString(wgsaLocationObject, "dbSNP"))
                    .cosmicId(optionalString(wgsaLocationObject, "COSMIC_ID"))
                    .phyloP46wayPlacental(string(wgsaLocationObject, "phyloP46way_placental"))
                    .phyloP100wayVertebrate(string(wgsaLocationObject, "phyloP100way_vertebrate"))
                    .siPhy29wayLogOdds(string(wgsaLocationObject, "SiPhy_29way_logOdds"))
                    .gwasSNP(optionalString(wgsaLocationObject, "GWAS_SNP"))
                    .gwasDIS(optionalString(wgsaLocationObject, "GWAS_DIS"))
                    .gwasPubmed(optionalString(wgsaLocationObject, "GWAS_PUBMED"))
                    .gerpRS(string(wgsaLocationObject, "GERP++_RS"))
                    .func(string(wgsaLocationObject, "Func"))
                    .wgRna(optionalString(wgsaLocationObject, "wgRna"))
                    .targetScanS(optionalString(wgsaLocationObject, "targetScanS"))
                    .key(string(wgsaLocationObject, "_key"))
                    .build());
        }
        return wgsaLocationList;
    }

    @NotNull
    private static List<MolecularMatchWGSAMap> createWGSAMaps(@Nullable JsonArray wgsaMapArray) {
        if (wgsaMapArray == null) {
            return Lists.newArrayList();
        }

        List<MolecularMatchWGSAMap> wgsaMapList = Lists.newArrayList();
        ViccDatamodelChecker wgsaMapChecker = ViccDatamodelCheckerFactory.molecularMatchWGSAMapChecker();

        for (JsonElement wgsaMapElement : wgsaMapArray) {
            JsonObject wgsaMapObject = wgsaMapElement.getAsJsonObject();
            wgsaMapChecker.check(wgsaMapObject);

            wgsaMapList.add(ImmutableMolecularMatchWGSAMap.builder()
                    .name(string(wgsaMapObject, "name"))
                    .gene(string(wgsaMapObject, "Gene"))
                    .transcript(string(wgsaMapObject, "Transcript"))
                    .exon(optionalString(wgsaMapObject, "Exon"))
                    .grch37ChrStartRefAlt(string(wgsaMapObject, "GRCh37_Chr_Start_Ref_Alt"))
                    .nucleotideChange(string(wgsaMapObject, "NucleotideChange"))
                    .aa(optionalString(wgsaMapObject, "AA"))
                    .synonyms(stringList(wgsaMapObject, "Synonyms"))
                    .protCoords(stringList(wgsaMapObject, "ProtCoords"))
                    .build());
        }
        return wgsaMapList;
    }

    @Nullable
    private static MolecularMatchExonsInfo createExonsInfo(@Nullable JsonObject exonsInfoObject) {
        if (exonsInfoObject == null) {
            return null;
        }

        ViccDatamodelCheckerFactory.molecularMatchExonsInfoChecker().check(exonsInfoObject);

        return ImmutableMolecularMatchExonsInfo.builder()
                .chr(string(exonsInfoObject, "chr"))
                .transcript(string(exonsInfoObject, "transcript"))
                .txStart(optionalString(exonsInfoObject, "txStart"))
                .txEnd(optionalString(exonsInfoObject, "txEnd"))
                .cdsStart(optionalString(exonsInfoObject, "cdsStart"))
                .cdsEnd(optionalString(exonsInfoObject, "cdsEnd"))
                .exonBoundaries(createExonBoundaries(exonsInfoObject.getAsJsonObject("exonBoundaries")))
                .build();
    }

    @NotNull
    private static MolecularMatchExonBoundaries createExonBoundaries(@NotNull JsonObject exonBoundariesObject) {
        ViccDatamodelCheckerFactory.molecularMatchExonBoundariesChecker().check(exonBoundariesObject);

        return ImmutableMolecularMatchExonBoundaries.builder()
                .exon1(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "1")))
                .exon2(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "2")))
                .exon3(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "3")))
                .exon4(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "4")))
                .exon5(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "5")))
                .exon6(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "6")))
                .exon7(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "7")))
                .exon8(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "8")))
                .exon9(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "9")))
                .exon10(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "10")))
                .exon11(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "11")))
                .exon12(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "12")))
                .exon13(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "13")))
                .exon14(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "14")))
                .exon15(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "15")))
                .exon16(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "16")))
                .exon17(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "17")))
                .exon18(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "18")))
                .exon19(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "19")))
                .exon20(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "20")))
                .exon21(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "21")))
                .exon22(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "22")))
                .exon23(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "23")))
                .exon24(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "24")))
                .exon25(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "25")))
                .exon26(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "26")))
                .exon27(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "27")))
                .exon28(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "28")))
                .exon29(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "29")))
                .exon30(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "30")))
                .exon31(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "31")))
                .exon32(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "32")))
                .exon33(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "33")))
                .exon34(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "34")))
                .exon35(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "35")))
                .exon36(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "36")))
                .exon37(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "37")))
                .exon38(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "38")))
                .exon39(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "39")))
                .exon40(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "40")))
                .exon41(createMolecularPosition(optionalJsonObject(exonBoundariesObject, "41")))
                .build();
    }

    @Nullable
    private static MolecularMatchPosition createMolecularPosition(@Nullable JsonObject positionObject) {
        if (positionObject == null) {
            return null;
        }

        ViccDatamodelCheckerFactory.molecularMatchPositionChecker().check(positionObject);

        return ImmutableMolecularMatchPosition.builder()
                .start(string(positionObject, "start"))
                .stop(string(positionObject, "stop"))
                .build();
    }

    @NotNull
    private static List<MolecularMatchFusionData> createFusionData(@Nullable JsonArray fusionDataArray) {
        if (fusionDataArray == null) {
            return Lists.newArrayList();
        }

        List<MolecularMatchFusionData> fusionDataList = Lists.newArrayList();
        ViccDatamodelChecker fusionDataChecker = ViccDatamodelCheckerFactory.molecularMatchFusionDataChecker();

        for (JsonElement fusionDataElement : fusionDataArray) {
            JsonObject fusionDataObject = fusionDataElement.getAsJsonObject();
            fusionDataChecker.check(fusionDataObject);

            fusionDataList.add(ImmutableMolecularMatchFusionData.builder()
                    .source(optionalString(fusionDataObject, "source"))
                    .synonym(optionalString(fusionDataObject, "synonym"))
                    .aChromosomes(optionalStringList(fusionDataObject, "Achr"))
                    .aBands(optionalStringList(fusionDataObject, "Aband"))
                    .aGenes(optionalStringList(fusionDataObject, "Agene"))
                    .aCoords(optionalStringList(fusionDataObject, "Acoord"))
                    .aTranscripts(optionalStringList(fusionDataObject, "Atx"))
                    .aOrientations(optionalStringList(fusionDataObject, "Aori"))
                    .aGenomicRegions(createFusionGenomicRegions(optionalJsonArray(fusionDataObject, "Agreg")))
                    .bChromosomes(optionalStringList(fusionDataObject, "Bchr"))
                    .bBands(optionalStringList(fusionDataObject, "Bband"))
                    .bGenes(optionalStringList(fusionDataObject, "Bgene"))
                    .bCoords(optionalStringList(fusionDataObject, "Bcoord"))
                    .bTranscripts(optionalStringList(fusionDataObject, "Btx"))
                    .bOrientations(optionalStringList(fusionDataObject, "Bori"))
                    .bGenomicRegions(createFusionGenomicRegions(optionalJsonArray(fusionDataObject, "Bgreg")))
                    .inserts(optionalStringList(fusionDataObject, "ins"))
                    .paper(optionalString(fusionDataObject, "Paper"))
                    .build());
        }
        return fusionDataList;
    }

    @NotNull
    private static List<MolecularMatchFusionGenomicRegion> createFusionGenomicRegions(@Nullable JsonArray genomicRegionArray) {
        if (genomicRegionArray == null) {
            return Lists.newArrayList();
        }

        List<MolecularMatchFusionGenomicRegion> genomicRegionList = Lists.newArrayList();
        ViccDatamodelChecker fusionGenomicRegionChecker = ViccDatamodelCheckerFactory.molecularMatchFusionGenomicRegionChecker();

        for (JsonElement genomicRegionElement : genomicRegionArray) {
            JsonObject genomicRegionObject = genomicRegionElement.getAsJsonObject();
            fusionGenomicRegionChecker.check(genomicRegionObject);

            genomicRegionList.add(ImmutableMolecularMatchFusionGenomicRegion.builder()
                    .num(string(genomicRegionObject, "num"))
                    .type(string(genomicRegionObject, "type"))
                    .build());
        }
        return genomicRegionList;
    }

    @NotNull
    private static List<MolecularMatchGRCh37Location> createGRCh37Locations(@NotNull JsonArray locationArray) {
        List<MolecularMatchGRCh37Location> grch37LocationList = Lists.newArrayList();
        ViccDatamodelChecker grch37LocationChecker = ViccDatamodelCheckerFactory.molecularMatchGRCh37LocationChecker();

        for (JsonElement locationElement : locationArray) {
            JsonObject locationObject = locationElement.getAsJsonObject();
            grch37LocationChecker.check(locationObject);

            grch37LocationList.add(ImmutableMolecularMatchGRCh37Location.builder()
                    .chr(nullableString(locationObject, "chr"))
                    .start(nullableString(locationObject, "start"))
                    .stop(nullableString(locationObject, "stop"))
                    .ref(nullableString(locationObject, "ref"))
                    .alt(nullableString(locationObject, "alt"))
                    .strand(string(locationObject, "strand"))
                    .transcriptConsequences(createGRCh37TranscriptConsequences(locationObject.getAsJsonArray("transcript_consequences")))
                    .validated(string(locationObject, "validated"))
                    .compositeKey(string(locationObject, "compositeKey"))
                    .build());
        }
        return grch37LocationList;
    }

    @NotNull
    private static List<MolecularMatchGRCh37TranscriptConsequence> createGRCh37TranscriptConsequences(
            @NotNull JsonArray transcriptConsequenceArray) {
        List<MolecularMatchGRCh37TranscriptConsequence> transcriptConsequenceList = Lists.newArrayList();
        ViccDatamodelChecker grch37TranscriptConsequenceChecker =
                ViccDatamodelCheckerFactory.molecularMatchGRCh37TranscriptConsequenceChecker();

        for (JsonElement transcriptConsequenceElement : transcriptConsequenceArray) {
            JsonObject transcriptConsequenceObject = transcriptConsequenceElement.getAsJsonObject();
            grch37TranscriptConsequenceChecker.check(transcriptConsequenceObject);

            transcriptConsequenceList.add(ImmutableMolecularMatchGRCh37TranscriptConsequence.builder()
                    .transcript(string(transcriptConsequenceObject, "transcript"))
                    .cdna(nullableString(transcriptConsequenceObject, "cdna"))
                    .aminoAcidChange(nullableString(transcriptConsequenceObject, "amino_acid_change"))
                    .txSites(stringList(transcriptConsequenceObject, "txSites"))
                    .intronNumber(nullableString(transcriptConsequenceObject, "intronNumber"))
                    .exonNumbers(stringList(transcriptConsequenceObject, "exonNumber"))
                    .build());
        }
        return transcriptConsequenceList;
    }

    @NotNull
    private static List<MolecularMatchVariantInfo> createVariantInfos(@NotNull JsonArray variantInfoArray) {
        List<MolecularMatchVariantInfo> variantInfoList = Lists.newArrayList();
        ViccDatamodelChecker variantInfoChecker = ViccDatamodelCheckerFactory.molecularMatchVariantInfoChecker();

        for (JsonElement variantInfoElement : variantInfoArray) {
            JsonObject variantInfoObject = variantInfoElement.getAsJsonObject();
            variantInfoChecker.check(variantInfoObject);

            variantInfoList.add(ImmutableMolecularMatchVariantInfo.builder()
                    .name(string(variantInfoObject, "name"))
                    .gene(string(variantInfoObject, "gene"))
                    .transcript(string(variantInfoObject, "transcript"))
                    .classification(string(variantInfoObject, "classification"))
                    .consequences(stringList(variantInfoObject, "consequences"))
                    .fusions(createFusions(variantInfoObject.getAsJsonArray("fusions")))
                    .locations(createLocations(variantInfoObject.getAsJsonArray("locations")))
                    .geneFusionPartner(string(variantInfoObject, "geneFusionPartner"))
                    .cosmicId(nullableString(variantInfoObject, "COSMIC_ID"))
                    .popFreqMax(string(variantInfoObject, "popFreqMax"))
                    .build());
        }
        return variantInfoList;
    }

    @NotNull
    private static List<MolecularMatchFusion> createFusions(@NotNull JsonArray fusionArray) {
        List<MolecularMatchFusion> fusionList = Lists.newArrayList();
        ViccDatamodelChecker fusionChecker = ViccDatamodelCheckerFactory.molecularMatchFusionChecker();

        for (JsonElement fusionElement : fusionArray) {
            JsonObject fusionObject = fusionElement.getAsJsonObject();
            fusionChecker.check(fusionObject);

            fusionList.add(ImmutableMolecularMatchFusion.builder()
                    .chr(string(fusionObject, "chr"))
                    .referenceGenome(string(fusionObject, "referenceGenome"))
                    .LBPWREP(string(fusionObject, "LBPWREP"))
                    .LBPWLEP(string(fusionObject, "LBPWLEP"))
                    .RBPWREP(string(fusionObject, "RBPWREP"))
                    .RBPWLEP(string(fusionObject, "RBPWLEP"))
                    .intronNumber(string(fusionObject, "intronNumber"))
                    .exonNumber(string(fusionObject, "exonNumber"))
                    .build());
        }
        return fusionList;
    }

    @NotNull
    private static List<MolecularMatchLocation> createLocations(@NotNull JsonArray locationArray) {
        List<MolecularMatchLocation> locationsList = Lists.newArrayList();
        ViccDatamodelChecker locationChecker = ViccDatamodelCheckerFactory.molecularMatchLocationChecker();

        for (JsonElement locationElement : locationArray) {
            JsonObject locationObject = locationElement.getAsJsonObject();
            locationChecker.check(locationObject);

            locationsList.add(ImmutableMolecularMatchLocation.builder()
                    .chr(string(locationObject, "chr"))
                    .start(string(locationObject, "start"))
                    .stop(string(locationObject, "stop"))
                    .ref(optionalString(locationObject, "ref"))
                    .alt(optionalString(locationObject, "alt"))
                    .cdna(optionalString(locationObject, "cdna"))
                    .aminoAcidChange(optionalString(locationObject, "amino_acid_change"))
                    .referenceGenome(optionalString(locationObject, "referenceGenome"))
                    .strand(optionalString(locationObject, "strand"))
                    .intronNumber(optionalString(locationObject, "intronNumber"))
                    .exonNumbers(optionalStringList(locationObject, "exonNumber"))
                    .build());
        }
        return locationsList;
    }

    @NotNull
    private static List<MolecularMatchPrevalence> createPrevalences(@NotNull JsonArray prevalenceArray) {
        List<MolecularMatchPrevalence> prevalenceList = Lists.newArrayList();
        ViccDatamodelChecker prevalenceChecker = ViccDatamodelCheckerFactory.molecularMatchPrevalenceChecker();

        for (JsonElement prevalenceElement : prevalenceArray) {
            JsonObject prevalenceObject = prevalenceElement.getAsJsonObject();
            prevalenceChecker.check(prevalenceObject);

            prevalenceList.add(ImmutableMolecularMatchPrevalence.builder()
                    .studyId(string(prevalenceObject, "studyId"))
                    .count(string(prevalenceObject, "count"))
                    .samples(string(prevalenceObject, "samples"))
                    .percent(string(prevalenceObject, "percent"))
                    .molecular(optionalString(prevalenceObject, "molecular"))
                    .condition(optionalString(prevalenceObject, "condition"))
                    .build());
        }
        return prevalenceList;
    }

    @NotNull
    private static List<MolecularMatchSource> createSources(@NotNull JsonArray sourceArray) {
        List<MolecularMatchSource> sourceList = Lists.newArrayList();
        ViccDatamodelChecker sourceChecker = ViccDatamodelCheckerFactory.molecularMatchSourceChecker();

        for (JsonElement sourceElement : sourceArray) {
            JsonObject sourceObject = sourceElement.getAsJsonObject();
            sourceChecker.check(sourceObject);

            sourceList.add(ImmutableMolecularMatchSource.builder()
                    .name(string(sourceObject, "name"))
                    .type(string(sourceObject, "type"))
                    .subType(optionalString(sourceObject, "subType"))
                    .valid(string(sourceObject, "valid"))
                    .pubId(string(sourceObject, "pubId"))
                    .link(string(sourceObject, "link"))
                    .trialId(optionalString(sourceObject, "trialId"))
                    .trialPhase(optionalString(sourceObject, "trialPhase"))
                    .year(string(sourceObject, "year"))
                    .functionalConsequence(optionalString(sourceObject, "functionalConsequence"))
                    .institution(optionalString(sourceObject, "institution"))
                    .trustRating(optionalString(sourceObject, "trustRating"))
                    .suppress(string(sourceObject, "suppress"))
                    .id(string(sourceObject, "id"))
                    .build());
        }
        return sourceList;
    }

    @NotNull
    private static List<MolecularMatchTierExplanation> createTierExplanations(@NotNull JsonArray tierExplanationArray) {
        List<MolecularMatchTierExplanation> tierExplanationList = Lists.newArrayList();
        ViccDatamodelChecker tierExplanationChecker = ViccDatamodelCheckerFactory.molecularMatchTierExplanationChecker();

        for (JsonElement tierExplanationElement : tierExplanationArray) {
            JsonObject tierExplanationObject = tierExplanationElement.getAsJsonObject();
            tierExplanationChecker.check(tierExplanationObject);

            tierExplanationList.add(ImmutableMolecularMatchTierExplanation.builder()
                    .tier(string(tierExplanationObject, "tier"))
                    .step(string(tierExplanationObject, "step"))
                    .message(string(tierExplanationObject, "message"))
                    .success(string(tierExplanationObject, "success"))
                    .build());
        }
        return tierExplanationList;
    }

    @NotNull
    private static List<MolecularMatchCriteriaUnmet> createCriteriaUnmets(@NotNull JsonArray criteriaUnmetArray) {
        List<MolecularMatchCriteriaUnmet> criteriaUnmetList = Lists.newArrayList();
        ViccDatamodelChecker criteriaUnmetChecker = ViccDatamodelCheckerFactory.molecularMatchCriteriaUnmetChecker();

        for (JsonElement criteriaUnmetElement : criteriaUnmetArray) {
            JsonObject criteriaUnmetObject = criteriaUnmetElement.getAsJsonObject();
            criteriaUnmetChecker.check(criteriaUnmetObject);

            criteriaUnmetList.add(ImmutableMolecularMatchCriteriaUnmet.builder()
                    .term(string(criteriaUnmetObject, "term"))
                    .filterType(string(criteriaUnmetObject, "filterType"))
                    .priority(string(criteriaUnmetObject, "priority"))
                    .facet(string(criteriaUnmetObject, "facet"))
                    .valid(optionalString(criteriaUnmetObject, "valid"))
                    .transcript(optionalString(criteriaUnmetObject, "transcript"))
                    .isNew(optionalString(criteriaUnmetObject, "isNew"))
                    .generatedBy(optionalString(criteriaUnmetObject, "generatedBy"))
                    .generatedByTerm(optionalString(criteriaUnmetObject, "generatedByTerm"))
                    .suppress(string(criteriaUnmetObject, "suppress"))
                    .manualSuppress(optionalString(criteriaUnmetObject, "manualSuppress"))
                    .primary(optionalString(criteriaUnmetObject, "primary"))
                    .compositeKey(string(criteriaUnmetObject, "compositeKey"))
                    .custom(optionalString(criteriaUnmetObject, "custom"))
                    .build());
        }
        return criteriaUnmetList;
    }

    @NotNull
    private static MolecularMatchAst createAst(@NotNull JsonObject astObject) {
        ViccDatamodelCheckerFactory.molecularMatchAstChecker().check(astObject);

        return ImmutableMolecularMatchAst.builder()
                .type(string(astObject, "type"))
                .raw(optionalString(astObject, "raw"))
                .value(optionalString(astObject, "value"))
                .operator(optionalString(astObject, "operator"))
                .left(createAstLeft(optionalJsonObject(astObject, "left")))
                .right(createAstRight(optionalJsonObject(astObject, "right")))
                .build();
    }

    @Nullable
    private static MolecularMatchAstLeft createAstLeft(@Nullable JsonObject astLeftObject) {
        if (astLeftObject == null) {
            return null;
        }

        ViccDatamodelCheckerFactory.molecularMatchAstChecker().check(astLeftObject);

        return ImmutableMolecularMatchAstLeft.builder()
                .type(string(astLeftObject, "type"))
                .raw(optionalString(astLeftObject, "raw"))
                .value(optionalString(astLeftObject, "value"))
                .operator(optionalString(astLeftObject, "operator"))
                .left(createAstLeftLeft(optionalJsonObject(astLeftObject, "left")))
                .right(createAstLeftRight(optionalJsonObject(astLeftObject, "right")))
                .build();
    }

    @Nullable
    private static MolecularMatchAstLeftLeft createAstLeftLeft(@Nullable JsonObject astLeftLeftObject) {
        if (astLeftLeftObject == null) {
            return null;
        }

        ViccDatamodelCheckerFactory.molecularMatchAstChecker().check(astLeftLeftObject);

        // We ignore deeper down "left + right". This continues recursively. If operator is present, there is a left/right that is ignored.
        return ImmutableMolecularMatchAstLeftLeft.builder()
                .type(string(astLeftLeftObject, "type"))
                .raw(optionalString(astLeftLeftObject, "raw"))
                .value(optionalString(astLeftLeftObject, "value"))
                .operator(optionalString(astLeftLeftObject, "operator"))
                .build();
    }

    @Nullable
    private static MolecularMatchAstLeftRight createAstLeftRight(@Nullable JsonObject astLeftRightObject) {
        if (astLeftRightObject == null) {
            return null;
        }

        ViccDatamodelCheckerFactory.molecularMatchAstChecker().check(astLeftRightObject);

        // We ignore deeper down "left + right". This continues recursively. If operator is present, there is a left/right that is ignored.
        return ImmutableMolecularMatchAstLeftRight.builder()
                .type(string(astLeftRightObject, "type"))
                .raw(optionalString(astLeftRightObject, "raw"))
                .value(optionalString(astLeftRightObject, "value"))
                .operator(optionalString(astLeftRightObject, "operator"))
                .build();
    }

    @Nullable
    private static MolecularMatchAstRight createAstRight(@Nullable JsonObject astRightObject) {
        if (astRightObject == null) {
            return null;
        }

        ViccDatamodelCheckerFactory.molecularMatchAstChecker().check(astRightObject);

        return ImmutableMolecularMatchAstRight.builder()
                .type(string(astRightObject, "type"))
                .raw(optionalString(astRightObject, "raw"))
                .value(optionalString(astRightObject, "value"))
                .operator(optionalString(astRightObject, "operator"))
                .left(createAstRightLeft(optionalJsonObject(astRightObject, "left")))
                .right(createAstRightRight(optionalJsonObject(astRightObject, "right")))
                .build();
    }

    @Nullable
    private static MolecularMatchAstRightLeft createAstRightLeft(@Nullable JsonObject astRightLeftObject) {
        if (astRightLeftObject == null) {
            return null;
        }

        ViccDatamodelCheckerFactory.molecularMatchAstChecker().check(astRightLeftObject);

        // We ignore deeper down "left + right". This continues recursively. If operator is present, there is a left/right that is ignored.
        return ImmutableMolecularMatchAstRightLeft.builder()
                .type(string(astRightLeftObject, "type"))
                .raw(optionalString(astRightLeftObject, "raw"))
                .value(optionalString(astRightLeftObject, "value"))
                .operator(optionalString(astRightLeftObject, "operator"))
                .build();
    }

    @Nullable
    private static MolecularMatchAstRightRight createAstRightRight(@Nullable JsonObject astRightRightObject) {
        if (astRightRightObject == null) {
            return null;
        }

        ViccDatamodelCheckerFactory.molecularMatchAstChecker().check(astRightRightObject);

        // We ignore deeper down "left + right". This continues recursively. If operator is present, there is a left/right that is ignored.
        return ImmutableMolecularMatchAstRightRight.builder()
                .type(string(astRightRightObject, "type"))
                .raw(optionalString(astRightRightObject, "raw"))
                .value(optionalString(astRightRightObject, "value"))
                .operator(optionalString(astRightRightObject, "operator"))
                .build();
    }

    @NotNull
    private static List<MolecularMatchTag> createTags(@NotNull JsonArray tagArray) {
        List<MolecularMatchTag> tagList = Lists.newArrayList();
        ViccDatamodelChecker tagChecker = ViccDatamodelCheckerFactory.molecularMatchTagChecker();

        for (JsonElement tagElement : tagArray) {
            JsonObject tagObject = tagElement.getAsJsonObject();
            tagChecker.check(tagObject);

            tagList.add(ImmutableMolecularMatchTag.builder()
                    .term(string(tagObject, "term"))
                    .facet(string(tagObject, "facet"))
                    .filterType(optionalString(tagObject, "filterType"))
                    .priority(string(tagObject, "priority"))
                    .transcript(optionalString(tagObject, "transcript"))
                    .valid(optionalString(tagObject, "valid"))
                    .generatedBy(optionalString(tagObject, "generatedBy"))
                    .generatedByTerm(optionalString(tagObject, "generatedByTerm"))
                    .isNew(optionalString(tagObject, "isNew"))
                    .primary(optionalString(tagObject, "primary"))
                    .custom(optionalString(tagObject, "custom"))
                    .suppress(optionalString(tagObject, "suppress"))
                    .manualSuppress(optionalString(tagObject, "manualSuppress"))
                    .composite(optionalString(tagObject, "composite"))
                    .compositeKey(optionalString(tagObject, "compositeKey"))
                    .build());
        }
        return tagList;
    }

    @NotNull
    private static List<MolecularMatchClassification> createClassifications(@NotNull JsonArray classificationArray) {
        List<MolecularMatchClassification> classificationList = Lists.newArrayList();
        ViccDatamodelChecker classificationChecker = ViccDatamodelCheckerFactory.molecularMatchClassificationChecker();

        for (JsonElement classificationElement : classificationArray) {
            JsonObject classificationObject = classificationElement.getAsJsonObject();
            classificationChecker.check(classificationObject);

            classificationList.add(ImmutableMolecularMatchClassification.builder()
                    .name(optionalString(classificationObject, "name"))
                    .geneSymbol(optionalString(classificationObject, "geneSymbol"))
                    .expandGeneSearch(optionalString(classificationObject, "expandGeneSearch"))
                    .transcript(optionalNullableString(classificationObject, "transcript"))
                    .transcripts(optionalStringList(classificationObject, "transcripts"))
                    .chromosomes(optionalStringList(classificationObject, "Chr"))
                    .starts(optionalStringList(classificationObject, "Start"))
                    .ends(optionalStringList(classificationObject, "End"))
                    .refs(optionalStringList(classificationObject, "Ref"))
                    .alts(optionalStringList(classificationObject, "Alt"))
                    .nucleotideChanges(optionalStringList(classificationObject, "NucleotideChange"))
                    .exons(optionalStringList(classificationObject, "Exon"))
                    .exonicFuncs(optionalStringList(classificationObject, "ExonicFunc"))
                    .classification(string(classificationObject, "classification"))
                    .classificationOverride(optionalNullableString(classificationObject, "classificationOverride"))
                    .pathology(optionalStringList(classificationObject, "pathology"))
                    .copyNumberType(optionalNullableString(classificationObject, "copyNumberType"))
                    .drugsApprovedOnLabelCount(optionalString(classificationObject, "drugsApprovedOnLabelCount"))
                    .drugsApprovedOffLabelCount(optionalString(classificationObject, "drugsApprovedOffLabelCount"))
                    .drugsExperimentalCount(optionalString(classificationObject, "drugsExperimentalCount"))
                    .trialCount(optionalString(classificationObject, "trialCount"))
                    .publicationCount(optionalString(classificationObject, "publicationCount"))
                    .sources(optionalStringList(classificationObject, "sources"))
                    .dbSNPs(optionalStringList(classificationObject, "dbSNP"))
                    .cosmicIds(optionalStringList(classificationObject, "COSMIC_ID"))
                    .popFreqMaxes(optionalStringList(classificationObject, "PopFreqMax"))
                    .parents(createParents(optionalJsonArray(classificationObject, "parents")))
                    .rootTerm(optionalString(classificationObject, "rootTerm"))
                    .alias(optionalString(classificationObject, "alias"))
                    .priority(optionalString(classificationObject, "priority"))
                    .description(optionalString(classificationObject, "description"))
                    .build());
        }
        return classificationList;
    }

    @NotNull
    private static List<MolecularMatchParent> createParents(@Nullable JsonArray parentArray) {
        if (parentArray == null) {
            return Lists.newArrayList();
        }

        List<MolecularMatchParent> parentList = Lists.newArrayList();
        ViccDatamodelChecker parentChecker = ViccDatamodelCheckerFactory.molecularMatchParentChecker();

        for (JsonElement parentElement : parentArray) {
            JsonObject parentObject = parentElement.getAsJsonObject();
            parentChecker.check(parentObject);

            parentList.add(ImmutableMolecularMatchParent.builder()
                    .name(string(parentObject, "name"))
                    .type(optionalNullableString(parentObject, "type"))
                    .actionableParent(optionalString(parentObject, "actionableParent"))
                    .transcripts(stringList(parentObject, "transcripts"))
                    .build());
        }
        return parentList;
    }

    @NotNull
    private static List<MolecularMatchTherapeuticContext> createTherapeuticContexts(@NotNull JsonArray therapeuticContextArray) {
        List<MolecularMatchTherapeuticContext> therapeuticContextList = Lists.newArrayList();
        ViccDatamodelChecker therapeuticContextChecker = ViccDatamodelCheckerFactory.molecularMatchTherapeuticContextChecker();

        for (JsonElement therapeuticContextElement : therapeuticContextArray) {
            JsonObject therapeuticContextObject = therapeuticContextElement.getAsJsonObject();
            therapeuticContextChecker.check(therapeuticContextObject);

            therapeuticContextList.add(ImmutableMolecularMatchTherapeuticContext.builder()
                    .name(string(therapeuticContextObject, "name"))
                    .facet(string(therapeuticContextObject, "facet"))
                    .suppress(string(therapeuticContextObject, "suppress"))
                    .valid(optionalString(therapeuticContextObject, "valid"))
                    .build());
        }
        return therapeuticContextList;
    }
}
