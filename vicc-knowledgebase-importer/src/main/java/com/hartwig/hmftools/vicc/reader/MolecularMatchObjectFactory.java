package com.hartwig.hmftools.vicc.reader;

import static com.hartwig.hmftools.vicc.reader.JsonFunctions.nullableString;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.optionalJsonArray;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.optionalJsonObject;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.optionalNullableString;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.optionalString;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.optionalStringList;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.string;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.stringList;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.toStringList;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatch;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAst;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAstLeft;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAstRight;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAstRightLeft;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAstRightLeftLeft;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAstRightLeftRight;
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
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAstRight;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAstRightLeft;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAstRightLeftLeft;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAstRightLeftRight;
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

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class MolecularMatchObjectFactory {

    private static final Logger LOGGER = LogManager.getLogger(MolecularMatchObjectFactory.class);

    private static final List<Integer> EXPECTED_MOLECULARMATCH_AST_SIZES = Lists.newArrayList(3, 4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_LEFT_SIZES = Lists.newArrayList(3, 4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_RIGHT_SIZES = Lists.newArrayList(3, 4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_RIGHT_LEFT_SIZES = Lists.newArrayList(3, 4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_RIGHT_RIGHT_SIZES = Lists.newArrayList(3);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_RIGHT_LEFT_RIGHT_SIZES = Lists.newArrayList(3);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_CLASSIFICATION_SIZES = Lists.newArrayList(3, 29, 30, 31);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TAGS_SIZES = Lists.newArrayList(3, 8, 9, 12, 13);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_THERAPEUTIC_CONTEXT_SIZES = Lists.newArrayList(3, 4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_PARENTS_SIZES = Lists.newArrayList(3, 4);

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
                .includeGene1(stringList(molecularMatchObject, "includeGene1"))
                .includeFinding1(stringList(molecularMatchObject, "includeFinding1"))
                .includeCondition1(stringList(molecularMatchObject, "includeCondition1"))
                .includeMutation1(optionalStringList(molecularMatchObject, "includeMutation1"))
                .includeDrug1(optionalStringList(molecularMatchObject, "includeDrug1"))
                .includeDrugClass1(optionalStringList(molecularMatchObject, "includeDrugclass1"))
                .includeResistance1(optionalStringList(molecularMatchObject, "includeResistance1"))
                .includeStage0(optionalStringList(molecularMatchObject, "includeStage0"))
                .includeGene0(optionalStringList(molecularMatchObject, "includeDrug0"))
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
        Set<String> keysAst = astObject.keySet();
        if (!EXPECTED_MOLECULARMATCH_AST_SIZES.contains(keysAst.size())) {
            LOGGER.warn("Found {} in molecular match ast rather than the expected {}", keysAst.size(), EXPECTED_MOLECULARMATCH_AST_SIZES);
            LOGGER.warn(keysAst);
        }
        return ImmutableMolecularMatchAst.builder()
                .raw(astObject.get("raw") == null ? null : astObject.getAsJsonPrimitive("raw").getAsString())
                .value(!astObject.has("value") ? null : astObject.getAsJsonPrimitive("value").getAsString())
                .operator(!astObject.has("operator") ? null : astObject.getAsJsonPrimitive("operator").getAsString())
                .right(!astObject.has("right") ? null : createRight(astObject.getAsJsonObject("right")))
                .type(astObject.getAsJsonPrimitive("type").getAsString())
                .left(!astObject.has("left") ? null : createLeft(astObject.getAsJsonObject("left")))
                .build();
    }

    @NotNull
    private static MolecularMatchAstLeft createLeft(@NotNull JsonObject objectLeft) {
        Set<String> keysLeft = objectLeft.keySet();
        if (!EXPECTED_MOLECULARMATCH_LEFT_SIZES.contains(keysLeft.size())) {
            LOGGER.warn("Found {} in molecular match ast left rather than the expected {}",
                    keysLeft.size(),
                    EXPECTED_MOLECULARMATCH_LEFT_SIZES);
            LOGGER.warn(keysLeft);
        }

        return ImmutableMolecularMatchAstLeft.builder()
                .operator(!objectLeft.has("operator") ? null : objectLeft.getAsJsonPrimitive("operator").getAsString())
                .raw(!objectLeft.has("raw") ? null : objectLeft.getAsJsonPrimitive("raw").getAsString())
                .type(objectLeft.getAsJsonPrimitive("type").getAsString())
                .value(!objectLeft.has("value") ? null : objectLeft.getAsJsonPrimitive("value").getAsString())
                .build();
    }

    @NotNull
    private static MolecularMatchAstRight createRight(@NotNull JsonObject objectRight) {
        Set<String> keysRight = objectRight.keySet();
        if (!EXPECTED_MOLECULARMATCH_RIGHT_SIZES.contains(keysRight.size())) {
            LOGGER.warn("Found {} in molecular match ast right rather than the expected {}",
                    keysRight.size(),
                    EXPECTED_MOLECULARMATCH_RIGHT_SIZES);
            LOGGER.warn(keysRight);
        }

        return ImmutableMolecularMatchAstRight.builder()
                .operator(!objectRight.has("operator") ? null : objectRight.getAsJsonPrimitive("operator").getAsString())
                .left(!objectRight.has("left") ? null : createRightLeft(objectRight.getAsJsonObject("left")))
                .right(!objectRight.has("right") ? null : createRightRight(objectRight.getAsJsonObject("right")))
                .raw(!objectRight.has("raw") ? null : objectRight.getAsJsonPrimitive("raw").getAsString())
                .type(objectRight.getAsJsonPrimitive("type").getAsString())
                .value(!objectRight.has("value") ? null : objectRight.getAsJsonPrimitive("value").getAsString())
                .build();
    }

    @NotNull
    private static MolecularMatchAstRightLeft createRightLeft(@NotNull JsonObject objectRight) {
        Set<String> keysRightLeft = objectRight.keySet();
        if (!EXPECTED_MOLECULARMATCH_RIGHT_LEFT_SIZES.contains(keysRightLeft.size())) {
            LOGGER.warn("Found {} in molecular match ast right left rather than the expected {}",
                    keysRightLeft.size(),
                    EXPECTED_MOLECULARMATCH_RIGHT_LEFT_SIZES);
            LOGGER.warn(keysRightLeft);
        }

        return ImmutableMolecularMatchAstRightLeft.builder()
                .raw(!objectRight.has("raw") ? null : objectRight.getAsJsonPrimitive("raw").getAsString())
                .type(objectRight.getAsJsonPrimitive("type").getAsString())
                .value(!objectRight.has("value") ? null : objectRight.getAsJsonPrimitive("value").getAsString())
                .right(!objectRight.has("right") ? null : createRightLeftRight(objectRight.getAsJsonObject("right")))
                .left(!objectRight.has("left") ? null : createRightLeftLeft(objectRight.getAsJsonObject("left")))
                .build();
    }

    @NotNull
    private static MolecularMatchAstRightRight createRightRight(@NotNull JsonObject objectRight) {
        Set<String> keysRightRight = objectRight.keySet();
        if (!EXPECTED_MOLECULARMATCH_RIGHT_RIGHT_SIZES.contains(keysRightRight.size())) {
            LOGGER.warn("Found {} in molecular match ast right right rather than the expected {}",
                    keysRightRight.size(),
                    EXPECTED_MOLECULARMATCH_RIGHT_RIGHT_SIZES);
            LOGGER.warn(keysRightRight);
        }

        return ImmutableMolecularMatchAstRightRight.builder()
                .raw(!objectRight.has("raw") ? null : objectRight.getAsJsonPrimitive("raw").getAsString())
                .type(objectRight.getAsJsonPrimitive("type").getAsString())
                .value(!objectRight.has("value") ? null : objectRight.getAsJsonPrimitive("value").getAsString())
                .build();
    }

    @NotNull
    private static MolecularMatchAstRightLeftLeft createRightLeftLeft(@NotNull JsonObject objectLeft) {
        Set<String> keysLeftRight = objectLeft.keySet();
        if (!EXPECTED_MOLECULARMATCH_RIGHT_LEFT_RIGHT_SIZES.contains(keysLeftRight.size())) {
            LOGGER.warn("Found {} in molecular match ast right left right rather than the expected {}",
                    keysLeftRight.size(),
                    EXPECTED_MOLECULARMATCH_RIGHT_LEFT_RIGHT_SIZES);
            LOGGER.warn(keysLeftRight);
        }

        return ImmutableMolecularMatchAstRightLeftLeft.builder()
                .raw(!objectLeft.has("raw") ? null : objectLeft.getAsJsonPrimitive("raw").getAsString())
                .type(objectLeft.getAsJsonPrimitive("type").getAsString())
                .value(!objectLeft.has("value") ? null : objectLeft.getAsJsonPrimitive("value").getAsString())
                .build();
    }

    @NotNull
    private static MolecularMatchAstRightLeftRight createRightLeftRight(@NotNull JsonObject objectRight) {
        Set<String> keysLeftRight = objectRight.keySet();
        if (!EXPECTED_MOLECULARMATCH_RIGHT_LEFT_RIGHT_SIZES.contains(keysLeftRight.size())) {
            LOGGER.warn("Found {} in molecular match ast right left right rather than the expected {}",
                    keysLeftRight.size(),
                    EXPECTED_MOLECULARMATCH_RIGHT_LEFT_RIGHT_SIZES);
            LOGGER.warn(keysLeftRight);
        }

        return ImmutableMolecularMatchAstRightLeftRight.builder()
                .raw(!objectRight.has("raw") ? null : objectRight.getAsJsonPrimitive("raw").getAsString())
                .type(objectRight.getAsJsonPrimitive("type").getAsString())
                .value(!objectRight.has("value") ? null : objectRight.getAsJsonPrimitive("value").getAsString())
                .build();
    }

    @NotNull
    private static List<MolecularMatchTherapeuticContext> createTherapeuticContexts(@NotNull JsonArray arrayTherapeuticContext) {
        List<MolecularMatchTherapeuticContext> therapeuticContextList = Lists.newArrayList();
        for (JsonElement therapeuticContext : arrayTherapeuticContext) {
            Set<String> keysTherapeuticContext = therapeuticContext.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_THERAPEUTIC_CONTEXT_SIZES.contains(keysTherapeuticContext.size())) {
                LOGGER.warn("Found {} in molecular match therapeutic context rather than the expected {}",
                        keysTherapeuticContext.size(),
                        EXPECTED_MOLECULARMATCH_THERAPEUTIC_CONTEXT_SIZES);
                LOGGER.warn(keysTherapeuticContext);
            }

            therapeuticContextList.add(ImmutableMolecularMatchTherapeuticContext.builder()
                    .facet(therapeuticContext.getAsJsonObject().getAsJsonPrimitive("facet").getAsString())
                    .suppress(therapeuticContext.getAsJsonObject().getAsJsonPrimitive("suppress").getAsString())
                    .valid(!therapeuticContext.getAsJsonObject().has("valid")
                            ? null
                            : therapeuticContext.getAsJsonObject().getAsJsonPrimitive("valid").getAsString())
                    .name(therapeuticContext.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .build());
        }
        return therapeuticContextList;
    }

    @NotNull
    private static List<MolecularMatchClassification> createClassifications(@NotNull JsonArray objectClassifications) {
        List<MolecularMatchClassification> classificationList = Lists.newArrayList();
        for (JsonElement classification : objectClassifications) {
            Set<String> keysClassification = classification.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_CLASSIFICATION_SIZES.contains(keysClassification.size())) {
                LOGGER.warn("Found {} in molecular match classification rather than the expected {}",
                        keysClassification.size(),
                        EXPECTED_MOLECULARMATCH_CLASSIFICATION_SIZES);
                LOGGER.warn(keysClassification);
            }

            classificationList.add(ImmutableMolecularMatchClassification.builder()
                    .end(!classification.getAsJsonObject().has("End")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("End")))
                    .classification(classification.getAsJsonObject().getAsJsonPrimitive("classification").getAsString())
                    .classificationOverride(
                            !classification.getAsJsonObject().has("classificationOverride") || classification.getAsJsonObject()
                                    .get("classificationOverride")
                                    .isJsonNull()
                                    ? null
                                    : classification.getAsJsonObject().getAsJsonPrimitive("classificationOverride").getAsString())
                    .start(!classification.getAsJsonObject().has("Start")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("Start")))
                    .chr(!classification.getAsJsonObject().has("Chr")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("Chr")))
                    .geneSymbol(!classification.getAsJsonObject().has("geneSymbol")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("geneSymbol").getAsString())
                    .pathology(!classification.getAsJsonObject().has("pathology")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("pathology")))
                    .ref(!classification.getAsJsonObject().has("Ref")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("Ref")))
                    .description(!classification.getAsJsonObject().has("description")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("description").getAsString())
                    .priority(!classification.getAsJsonObject().has("priority")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("priority").getAsString())
                    .nucleotideChange(!classification.getAsJsonObject().has("NucleotideChange")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("NucleotideChange")))
                    .parents(!classification.getAsJsonObject().has("parents")
                            ? null
                            : createParents(classification.getAsJsonObject().getAsJsonArray("parents")))
                    .expandGeneSearch(!classification.getAsJsonObject().has("expandGeneSearch")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("expandGeneSearch").getAsString())
                    .drugsExperimentalCount(!classification.getAsJsonObject().has("drugsExperimentalCount")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("drugsExperimentalCount").getAsString())
                    .exon(!classification.getAsJsonObject().has("Exon")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("Exon")))
                    .drugsApprovedOffLabelCount(!classification.getAsJsonObject().has("drugsApprovedOffLabelCount")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("drugsApprovedOffLabelCount").getAsString())
                    .exonicFunc(!classification.getAsJsonObject().has("ExonicFunc")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("ExonicFunc")))
                    .popFreqMax(!classification.getAsJsonObject().has("PopFreqMax")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("PopFreqMax")))
                    .copyNumberType(!classification.getAsJsonObject().has("copyNumberType") || classification.getAsJsonObject()
                            .get("copyNumberType")
                            .isJsonNull() ? null : classification.getAsJsonObject().getAsJsonPrimitive("copyNumberType").getAsString())
                    .publicationCount(!classification.getAsJsonObject().has("publicationCount")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("publicationCount").getAsString())
                    .transcript(!classification.getAsJsonObject().has("transcript") || classification.getAsJsonObject()
                            .get("transcript")
                            .isJsonNull() ? null : classification.getAsJsonObject().getAsJsonPrimitive("transcript").getAsString())
                    .dbSNP(!classification.getAsJsonObject().has("dbSNP")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("dbSNP")))
                    .alt(!classification.getAsJsonObject().has("Alt")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("Alt")))
                    .name(!classification.getAsJsonObject().has("name")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .rootTerm(!classification.getAsJsonObject().has("rootTerm")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("rootTerm").getAsString())
                    .sources(!classification.getAsJsonObject().has("sources")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("sources")))
                    .drugsApprovedOnLabelCount(!classification.getAsJsonObject().has("drugsApprovedOnLabelCount")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("drugsApprovedOnLabelCount").getAsString())
                    .trialCount(!classification.getAsJsonObject().has("trialCount")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("trialCount").getAsString())
                    .alias(!classification.getAsJsonObject().has("alias")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("alias").getAsString())
                    .cosmicId(!classification.getAsJsonObject().has("COSMIC_ID")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("COSMIC_ID")))
                    .transcripts(!classification.getAsJsonObject().has("transcripts")
                            ? null
                            : toStringList(classification.getAsJsonObject().getAsJsonArray("transcripts")))
                    .build());
        }
        return classificationList;
    }

    @NotNull
    private static List<MolecularMatchParent> createParents(@NotNull JsonArray arrayParents) {
        List<MolecularMatchParent> parentsList = Lists.newArrayList();
        for (JsonElement parents : arrayParents) {
            Set<String> keysParents = parents.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_PARENTS_SIZES.contains(keysParents.size())) {
                LOGGER.warn("Found {} in molecular match parents rather than the expected {}",
                        keysParents.size(),
                        EXPECTED_MOLECULARMATCH_PARENTS_SIZES);
                LOGGER.warn(keysParents);
            }

            parentsList.add(ImmutableMolecularMatchParent.builder()
                    .transcripts(toStringList(parents.getAsJsonObject().getAsJsonArray("transcripts")))
                    .type(!parents.getAsJsonObject().has("type") || parents.getAsJsonObject().get("type").isJsonNull()
                            ? null
                            : parents.getAsJsonObject().getAsJsonPrimitive("type").getAsString())
                    .name(parents.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .actionableParent(!parents.getAsJsonObject().has("actionableParent")
                            ? null
                            : parents.getAsJsonObject().getAsJsonPrimitive("actionableParent").getAsString())
                    .build());
        }
        return parentsList;
    }

    @NotNull
    private static List<MolecularMatchTag> createTags(@NotNull JsonArray arrayTags) {
        List<MolecularMatchTag> tagsList = Lists.newArrayList();
        for (JsonElement tags : arrayTags) {
            Set<String> keysTags = tags.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_TAGS_SIZES.contains(keysTags.size())) {
                LOGGER.warn("Found {} in molecular match tags rather than the expected {}",
                        keysTags.size(),
                        EXPECTED_MOLECULARMATCH_TAGS_SIZES);
                LOGGER.warn(keysTags);
            }

            tagsList.add(ImmutableMolecularMatchTag.builder()
                    .priority(tags.getAsJsonObject().getAsJsonPrimitive("priority").getAsString())
                    .compositeKey(!tags.getAsJsonObject().has("compositeKey")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("compositeKey").getAsString())
                    .suppress(!tags.getAsJsonObject().has("suppress")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("suppress").getAsString())
                    .filterType(!tags.getAsJsonObject().has("filterType")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("filterType").getAsString())
                    .term(tags.getAsJsonObject().getAsJsonPrimitive("term").getAsString())
                    .primary(!tags.getAsJsonObject().has("primary")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("primary").getAsString())
                    .facet(tags.getAsJsonObject().getAsJsonPrimitive("facet").getAsString())
                    .valid(!tags.getAsJsonObject().has("valid") ? null : tags.getAsJsonObject().getAsJsonPrimitive("valid").getAsString())
                    .custom(!tags.getAsJsonObject().has("custom")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("custom").getAsString())
                    .isNew(!tags.getAsJsonObject().has("isNew") ? null : tags.getAsJsonObject().getAsJsonPrimitive("isNew").getAsString())
                    .generatedBy(!tags.getAsJsonObject().has("generatedBy")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("generatedBy").getAsString())
                    .manualSuppress(!tags.getAsJsonObject().has("manualSuppress")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("manualSuppress").getAsString())
                    .generatedByTerm(!tags.getAsJsonObject().has("generatedByTerm")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("generatedByTerm").getAsString())
                    .transcript(!tags.getAsJsonObject().has("transcript")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("transcript").getAsString())
                    .build());
        }
        return tagsList;
    }
}
