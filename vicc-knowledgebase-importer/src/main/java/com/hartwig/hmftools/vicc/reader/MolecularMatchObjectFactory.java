package com.hartwig.hmftools.vicc.reader;

import static com.hartwig.hmftools.vicc.reader.JsonFunctions.toStringList;

import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatch;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAreg;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAst;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAstLeft;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAstRight;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAstRightLeft;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAstRightLeftRight;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchAstRightRight;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchBreg;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchClassification;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchCriteriaUnmet;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchExonInfo;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchExonsBoundries;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchFusionData;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchFusions;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchGRch37Location;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchLocations;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchMutations;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchParents;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchPositions;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchPrevalence;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchSource;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchTags;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchTherapeuticContext;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchTierExplanation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchTranscriptConsequence;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchTranscriptConsequencesGRCH37;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchVariantInfo;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchWGSaMap;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.ImmutableMolecularMatchWGSadataLocation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatch;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAreg;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAst;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAstLeft;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAstRight;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAstRightLeft;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAstRightLeftRight;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchAstRightRight;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchBreg;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchClassification;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchCriteriaUnmet;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchExonInfo;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchExonsBoundries;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchFusionData;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchFusions;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchGRch37Location;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchLocations;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchMutations;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchParents;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchPositions;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchPrevalence;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchSource;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchTags;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchTherapeuticContext;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchTierExplanation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchTranscriptConsequence;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchTranscriptConsequencesGRCH37;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchVariantInfo;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchWGSaMap;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatchWGSadataLocation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class MolecularMatchObjectFactory {

    private static final Logger LOGGER = LogManager.getLogger(MolecularMatchObjectFactory.class);

    private static final List<Integer> EXPECTED_MOLECULARMATCH_ELEMENT_SIZES = Lists.newArrayList(34, 35, 36, 37, 38, 39, 40, 41, 42);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_AST_SIZES = Lists.newArrayList(3, 4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_LEFT_SIZES = Lists.newArrayList(3, 4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_RIGHT_SIZES = Lists.newArrayList(3, 4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_RIGHT_LEFT_SIZES = Lists.newArrayList(3, 4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_RIGHT_RIGHT_SIZES = Lists.newArrayList(3);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_RIGHT_LEFT_RIGHT_SIZES = Lists.newArrayList(3);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_CLASSIFICATION_SIZES = Lists.newArrayList(3, 29, 30, 31);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_CRITERIA_UNMET_SIZES = Lists.newArrayList(8, 9, 12, 13);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_LOCATIONGRCH37_SIZES = Lists.newArrayList(9);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_LOCATIONS_SIZES = Lists.newArrayList(3, 11);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_MUTATIONS_SIZES = Lists.newArrayList(13, 14, 16, 17, 18, 19);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_PREVALENCE_SIZES = Lists.newArrayList(4, 6);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_SOURCE_SIZES = Lists.newArrayList(8, 9, 10, 11, 12, 13);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TAGS_SIZES = Lists.newArrayList(3, 8, 9, 12, 13);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_THERAPEUTIC_CONTEXT_SIZES = Lists.newArrayList(3, 4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TIER_EXPLANATION_SIZES = Lists.newArrayList(4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TRANSCRIPT_CONSEQUENCES_SIZES = Lists.newArrayList(9, 14, 15, 16);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TRANSCRIPT_CONSEQUENCES__GRCH37_SIZES = Lists.newArrayList(6);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_VARIANTINFO_SIZES = Lists.newArrayList(10);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_PARENTS_SIZES = Lists.newArrayList(3, 4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_FUSIONS_SIZES = Lists.newArrayList(8);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_WGSADATA_SIZES = Lists.newArrayList(1);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_WGSAMAP_SIZES = Lists.newArrayList(7, 9);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_WGSADATA_LOCATION_SIZES =
            Lists.newArrayList(20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 39, 40, 41, 45);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_EXONSINFO_SIZES = Lists.newArrayList(3, 7);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_EXONSBOUNDARIES_SIZES =
            Lists.newArrayList(1, 2, 3, 5, 6, 7, 8, 9, 11, 13, 16, 17, 20, 21, 22, 24, 26, 27, 28, 29, 38, 41);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_POSITIONS_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_FUSIONDATA_SIZES = Lists.newArrayList(1, 2, 15, 17);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_BREG_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_AREG_SIZES = Lists.newArrayList(2);

    private MolecularMatchObjectFactory() {
    }

    @NotNull
    static MolecularMatch create(@NotNull JsonObject objectMolecularMatch) {
        Set<String> keysMolecularMatch = objectMolecularMatch.keySet();
        if (!EXPECTED_MOLECULARMATCH_ELEMENT_SIZES.contains(keysMolecularMatch.size())) {
            LOGGER.warn("Found {} in molecular match rather than the expected {}",
                    keysMolecularMatch.size(),
                    EXPECTED_MOLECULARMATCH_ELEMENT_SIZES);
            LOGGER.warn(keysMolecularMatch);
        }

        return ImmutableMolecularMatch.builder()
                .criteriaUnmet(createCriteriaUnmet(objectMolecularMatch.getAsJsonArray("criteriaUnmet")))
                .prevalence(createPrevalence(objectMolecularMatch.getAsJsonArray("prevalence")))
                .score(objectMolecularMatch.getAsJsonPrimitive("_score").getAsString())
                .autoGenerateNarrative(objectMolecularMatch.getAsJsonPrimitive("autoGenerateNarrative").getAsString())
                .mutations(createMutations(objectMolecularMatch.getAsJsonArray("mutations")))
                .sources(createSource(objectMolecularMatch.getAsJsonArray("sources")))
                .clinicalSignificance(objectMolecularMatch.getAsJsonPrimitive("clinicalSignificance").getAsString())
                .id(objectMolecularMatch.getAsJsonPrimitive("id").getAsString())
                .includeCondition0(toStringList(objectMolecularMatch.getAsJsonArray("includeCondition0")))
                .includeCondition1(toStringList(objectMolecularMatch.getAsJsonArray("includeCondition1")))
                .uniqueKey(objectMolecularMatch.getAsJsonPrimitive("uniqueKey").getAsString())
                .civicValue(objectMolecularMatch.getAsJsonPrimitive("civic").getAsString())
                .hashKey(objectMolecularMatch.getAsJsonPrimitive("hashKey").getAsString())
                .regulatoryBodyApproved(objectMolecularMatch.getAsJsonPrimitive("regulatoryBodyApproved").getAsString())
                .version(objectMolecularMatch.getAsJsonPrimitive("version").getAsString())
                .includeMutation1(!objectMolecularMatch.has("includeMutation1")
                        ? null
                        : toStringList(objectMolecularMatch.getAsJsonArray("includeMutation1")))
                .includeMutation0(!objectMolecularMatch.has("includeMutation0")
                        ? null
                        : toStringList(objectMolecularMatch.getAsJsonArray("includeMutation0")))
                .guidelineBody(!objectMolecularMatch.has("guidelineBody")
                        ? null
                        : objectMolecularMatch.getAsJsonPrimitive("guidelineBody").getAsString())
                .regulatoryBody(objectMolecularMatch.getAsJsonPrimitive("regulatoryBody").getAsString())
                .customer(objectMolecularMatch.getAsJsonPrimitive("customer").getAsString())
                .direction(objectMolecularMatch.getAsJsonPrimitive("direction").getAsString())
                .ampcap(objectMolecularMatch.getAsJsonPrimitive("ampcap").getAsString())
                .asts(createAst(objectMolecularMatch.getAsJsonObject("ast")))
                .variantInfo(createVariantInfo(objectMolecularMatch.getAsJsonArray("variantInfo")))
                .guidelineVersion(!objectMolecularMatch.has("guidelineVersion")
                        ? null
                        : objectMolecularMatch.getAsJsonPrimitive("guidelineVersion").getAsString())
                .institution(!objectMolecularMatch.has("institution")
                        ? null
                        : toStringList(objectMolecularMatch.getAsJsonArray("institution")))
                .tier(objectMolecularMatch.getAsJsonPrimitive("tier").getAsString())
                .tierExplanation(createTierExplanation(objectMolecularMatch.getAsJsonArray("tierExplanation")))
                .mvld(objectMolecularMatch.getAsJsonPrimitive("mvld").getAsString())
                .tags(createTags(objectMolecularMatch.getAsJsonArray("tags")))
                .criteriaMet(toStringList(objectMolecularMatch.getAsJsonArray("criteriaMet")))
                .biomarkerClass(objectMolecularMatch.getAsJsonPrimitive("biomarkerClass").getAsString())
                .classification(createClassification(objectMolecularMatch.getAsJsonArray("classifications")))
                .includeDrug1(!objectMolecularMatch.has("includeDrug1")
                        ? null
                        : toStringList(objectMolecularMatch.getAsJsonArray("includeDrug1")))
                .includeStage0(!objectMolecularMatch.has("includeStage0")
                        ? null
                        : toStringList(objectMolecularMatch.getAsJsonArray("includeStage0")))
                .therapeuticContext(createTherapeuticContext(objectMolecularMatch.getAsJsonArray("therapeuticContext")))
                .sixtier(objectMolecularMatch.getAsJsonPrimitive("sixtier").getAsString())
                .noTherapyAvailable(!objectMolecularMatch.has("noTherapyAvailable")
                        ? null
                        : objectMolecularMatch.getAsJsonPrimitive("noTherapyAvailable").getAsString())
                .external_id(!objectMolecularMatch.has("external_id")
                        ? null
                        : toStringList(objectMolecularMatch.getAsJsonArray("external_id")))
                .narrative(objectMolecularMatch.getAsJsonPrimitive("narrative").getAsString())
                .expression(objectMolecularMatch.getAsJsonPrimitive("expression").getAsString())
                .includeGene0(!objectMolecularMatch.has("includeDrug0")
                        ? null
                        : toStringList(objectMolecularMatch.getAsJsonArray("includeGene0")))
                .build();
    }

    @NotNull
    private static List<MolecularMatchTherapeuticContext> createTherapeuticContext(@NotNull JsonArray arrayTherapeuticContext) {
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
    private static List<MolecularMatchClassification> createClassification(@NotNull JsonArray objectClassifications) {
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
                    .NucleotideChange(!classification.getAsJsonObject().has("NucleotideChange")
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
                    .COSMIC_ID(!classification.getAsJsonObject().has("COSMIC_ID")
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
    private static List<MolecularMatchParents> createParents(@NotNull JsonArray arrayParents) {
        List<MolecularMatchParents> parentsList = Lists.newArrayList();
        for (JsonElement parents : arrayParents) {
            Set<String> keysParents = parents.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_PARENTS_SIZES.contains(keysParents.size())) {
                LOGGER.warn("Found {} in molecular match parents rather than the expected {}",
                        keysParents.size(),
                        EXPECTED_MOLECULARMATCH_PARENTS_SIZES);
                LOGGER.warn(keysParents);
            }

            parentsList.add(ImmutableMolecularMatchParents.builder()
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
    private static List<MolecularMatchTags> createTags(@NotNull JsonArray arrayTags) {
        List<MolecularMatchTags> tagsList = Lists.newArrayList();
        for (JsonElement tags : arrayTags) {
            Set<String> keysTags = tags.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_TAGS_SIZES.contains(keysTags.size())) {
                LOGGER.warn("Found {} in molecular match tags rather than the expected {}",
                        keysTags.size(),
                        EXPECTED_MOLECULARMATCH_TAGS_SIZES);
                LOGGER.warn(keysTags);
            }

            tagsList.add(ImmutableMolecularMatchTags.builder()
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

    @NotNull
    private static List<MolecularMatchTierExplanation> createTierExplanation(@NotNull JsonArray arrarTierExplanation) {
        List<MolecularMatchTierExplanation> tierExplanationList = Lists.newArrayList();
        for (JsonElement tierExplanation : arrarTierExplanation) {
            Set<String> keysTierExplanation = tierExplanation.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_TIER_EXPLANATION_SIZES.contains(keysTierExplanation.size())) {
                LOGGER.warn("Found {} in molecular match tier explanation rather than the expected {}",
                        keysTierExplanation.size(),
                        EXPECTED_MOLECULARMATCH_TIER_EXPLANATION_SIZES);
                LOGGER.warn(keysTierExplanation);
            }

            tierExplanationList.add(ImmutableMolecularMatchTierExplanation.builder()
                    .tier(tierExplanation.getAsJsonObject().getAsJsonPrimitive("tier").getAsString())
                    .step(tierExplanation.getAsJsonObject().getAsJsonPrimitive("step").getAsString())
                    .message(tierExplanation.getAsJsonObject().getAsJsonPrimitive("message").getAsString())
                    .success(tierExplanation.getAsJsonObject().getAsJsonPrimitive("success").getAsString())
                    .build());
        }
        return tierExplanationList;
    }

    @NotNull
    private static List<MolecularMatchVariantInfo> createVariantInfo(@NotNull JsonArray arrayVariantInfo) {
        List<MolecularMatchVariantInfo> variantInfoList = Lists.newArrayList();

        for (JsonElement variantInfo : arrayVariantInfo) {
            Set<String> keysVariantInfo = variantInfo.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_VARIANTINFO_SIZES.contains(keysVariantInfo.size())) {
                LOGGER.warn("Found {} in molecular match variant info rather than the expected {}",
                        keysVariantInfo.size(),
                        EXPECTED_MOLECULARMATCH_VARIANTINFO_SIZES);
                LOGGER.warn(keysVariantInfo);
            }

            variantInfoList.add(ImmutableMolecularMatchVariantInfo.builder()
                    .classification(variantInfo.getAsJsonObject().getAsJsonPrimitive("classification").getAsString())
                    .name(variantInfo.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .consequences(toStringList(variantInfo.getAsJsonObject().getAsJsonArray("consequences")))
                    .fusions(createFusions(variantInfo.getAsJsonObject().getAsJsonArray("fusions")))
                    .locations(createLocations(variantInfo.getAsJsonObject().getAsJsonArray("locations")))
                    .geneFusionPartner(variantInfo.getAsJsonObject().getAsJsonPrimitive("geneFusionPartner").getAsString())
                    .COSMIC_ID(variantInfo.getAsJsonObject().get("COSMIC_ID").isJsonNull()
                            ? null
                            : variantInfo.getAsJsonObject().getAsJsonPrimitive("COSMIC_ID").getAsString())
                    .gene(variantInfo.getAsJsonObject().getAsJsonPrimitive("gene").getAsString())
                    .transcript(variantInfo.getAsJsonObject().getAsJsonPrimitive("transcript").getAsString())
                    .popFreqMax(variantInfo.getAsJsonObject().getAsJsonPrimitive("popFreqMax").getAsString())
                    .build());
        }
        return variantInfoList;
    }

    @NotNull
    private static List<MolecularMatchFusions> createFusions(@NotNull JsonArray arrayFusions) {
        List<MolecularMatchFusions> fusionsList = Lists.newArrayList();

        for (JsonElement fusions : arrayFusions) {
            Set<String> keysFusions = fusions.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_FUSIONS_SIZES.contains(keysFusions.size())) {
                LOGGER.warn("Found {} in molecular match fusions rather than the expected {}",
                        keysFusions.size(),
                        EXPECTED_MOLECULARMATCH_FUSIONS_SIZES);
                LOGGER.warn(keysFusions);
            }
            fusionsList.add(ImmutableMolecularMatchFusions.builder()
                    .referenceGenome(fusions.getAsJsonObject().get("referenceGenome").getAsString())
                    .LBPWREP(fusions.getAsJsonObject().get("LBPWREP").getAsString())
                    .RBPWREP(fusions.getAsJsonObject().get("RBPWREP").getAsString())
                    .exonNumber(fusions.getAsJsonObject().get("exonNumber").getAsString())
                    .chr(fusions.getAsJsonObject().get("chr").getAsString())
                    .RBPWLEP(fusions.getAsJsonObject().get("RBPWLEP").getAsString())
                    .intronNumber(fusions.getAsJsonObject().get("intronNumber").getAsString())
                    .LBPWLEP(fusions.getAsJsonObject().get("LBPWLEP").getAsString())
                    .build());
        }
        return fusionsList;
    }

    @NotNull
    private static List<MolecularMatchLocations> createLocations(@NotNull JsonArray arrayLocations) {
        List<MolecularMatchLocations> locationsList = Lists.newArrayList();
        for (JsonElement locations : arrayLocations) {
            Set<String> keysLocations = locations.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_LOCATIONS_SIZES.contains(keysLocations.size())) {
                LOGGER.warn("Found {} in molecular match locations rather than the expected {}",
                        keysLocations.size(),
                        EXPECTED_MOLECULARMATCH_LOCATIONS_SIZES);
                LOGGER.warn(keysLocations);
            }

            locationsList.add(ImmutableMolecularMatchLocations.builder()
                    .aminoAcidChange(!locations.getAsJsonObject().has("amino_acid_change")
                            ? null
                            : locations.getAsJsonObject().getAsJsonPrimitive("amino_acid_change").getAsString())
                    .intronNumber(!locations.getAsJsonObject().has("intronNumber")
                            ? null
                            : locations.getAsJsonObject().getAsJsonPrimitive("intronNumber").getAsString())
                    .exonNumber(!locations.getAsJsonObject().has("exonNumber") ? null : createArrayExonNumber(locations))
                    .stop(locations.getAsJsonObject().getAsJsonPrimitive("stop").getAsString())
                    .start(locations.getAsJsonObject().getAsJsonPrimitive("start").getAsString())
                    .chr(locations.getAsJsonObject().getAsJsonPrimitive("chr").getAsString())
                    .strand(!locations.getAsJsonObject().has("strand")
                            ? null
                            : locations.getAsJsonObject().getAsJsonPrimitive("strand").getAsString())
                    .alt(!locations.getAsJsonObject().has("alt")
                            ? null
                            : locations.getAsJsonObject().getAsJsonPrimitive("alt").getAsString())
                    .referenceGenome(!locations.getAsJsonObject().has("referenceGenome")
                            ? null
                            : locations.getAsJsonObject().getAsJsonPrimitive("referenceGenome").getAsString())
                    .ref(!locations.getAsJsonObject().has("ref")
                            ? null
                            : locations.getAsJsonObject().getAsJsonPrimitive("ref").getAsString())
                    .cdna(!locations.getAsJsonObject().has("cdna")
                            ? null
                            : locations.getAsJsonObject().getAsJsonPrimitive("cdna").getAsString())
                    .build());
        }
        return locationsList;
    }

    @NotNull
    private static List<String> createArrayExonNumber(@NotNull JsonElement elementExonNumber) {
        if (elementExonNumber.getAsJsonObject().get("exonNumber").isJsonArray()) {
            return toStringList(elementExonNumber.getAsJsonObject().getAsJsonArray("exonNumber"));
        } else {
            return Collections.singletonList(elementExonNumber.getAsJsonObject().getAsJsonPrimitive("exonNumber").getAsString());
        }
    }

    @NotNull
    private static MolecularMatchAst createAst(@NotNull JsonObject objectAst) {
        Set<String> keysAst = objectAst.keySet();
        if (!EXPECTED_MOLECULARMATCH_AST_SIZES.contains(keysAst.size())) {
            LOGGER.warn("Found {} in molecular match ast rather than the expected {}", keysAst.size(), EXPECTED_MOLECULARMATCH_AST_SIZES);
            LOGGER.warn(keysAst);
        }
        return ImmutableMolecularMatchAst.builder()
                .raw(objectAst.get("raw") == null ? null : objectAst.getAsJsonPrimitive("raw").getAsString())
                .value(!objectAst.has("value") ? null : objectAst.getAsJsonPrimitive("value").getAsString())
                .operator(!objectAst.has("operator") ? null : objectAst.getAsJsonPrimitive("operator").getAsString())
                .right(!objectAst.has("right") ? null : createRight(objectAst.getAsJsonObject("right")))
                .type(objectAst.getAsJsonPrimitive("type").getAsString())
                .left(!objectAst.has("left") ? null : createLeft(objectAst.getAsJsonObject("left")))
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
                .left(!objectRight.has("left") ? null : createLeft(objectRight.getAsJsonObject("left")))
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
    private static List<MolecularMatchSource> createSource(@NotNull JsonArray arraySources) {
        List<MolecularMatchSource> sourcesList = Lists.newArrayList();
        for (JsonElement source : arraySources) {
            Set<String> keysSource = source.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_SOURCE_SIZES.contains(keysSource.size())) {
                LOGGER.warn("Found {} in molecular match source rather than the expected {}",
                        keysSource.size(),
                        EXPECTED_MOLECULARMATCH_SOURCE_SIZES);
                LOGGER.warn(keysSource);
            }

            sourcesList.add(ImmutableMolecularMatchSource.builder()
                    .name(source.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .suppress(source.getAsJsonObject().getAsJsonPrimitive("suppress").getAsString())
                    .pubId(source.getAsJsonObject().getAsJsonPrimitive("pubId").getAsString())
                    .subType(!source.getAsJsonObject().has("subType")
                            ? null
                            : source.getAsJsonObject().getAsJsonPrimitive("subType").getAsString())
                    .valid(source.getAsJsonObject().getAsJsonPrimitive("valid").getAsString())
                    .link(source.getAsJsonObject().getAsJsonPrimitive("link").getAsString())
                    .year(source.getAsJsonObject().getAsJsonPrimitive("year").getAsString())
                    .trialId(!source.getAsJsonObject().has("trialId")
                            ? null
                            : source.getAsJsonObject().getAsJsonPrimitive("trialId").getAsString())
                    .type(source.getAsJsonObject().getAsJsonPrimitive("type").getAsString())
                    .id(source.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .institution(!source.getAsJsonObject().has("institution")
                            ? null
                            : source.getAsJsonObject().getAsJsonPrimitive("institution").getAsString())
                    .trialPhase(!source.getAsJsonObject().has("trialPhase")
                            ? null
                            : source.getAsJsonObject().getAsJsonPrimitive("trialPhase").getAsString())
                    .functionalConsequence(!source.getAsJsonObject().has("functionalConsequence")
                            ? null
                            : source.getAsJsonObject().getAsJsonPrimitive("functionalConsequence").getAsString())
                    .trustRating(!source.getAsJsonObject().has("trustRating")
                            ? null
                            : source.getAsJsonObject().getAsJsonPrimitive("trustRating").getAsString())
                    .build());
        }
        return sourcesList;
    }

    @NotNull
    private static List<MolecularMatchCriteriaUnmet> createCriteriaUnmet(@NotNull JsonArray arrayCriteriaUnmet) {
        List<MolecularMatchCriteriaUnmet> criteriaUnmetList = Lists.newArrayList();

        for (JsonElement criteriaUnmet : arrayCriteriaUnmet) {
            Set<String> keysCriteriaUnmet = criteriaUnmet.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_CRITERIA_UNMET_SIZES.contains(keysCriteriaUnmet.size())) {
                LOGGER.warn("Found {} in molecular match criteria unmet rather than the expected {}",
                        keysCriteriaUnmet.size(),
                        EXPECTED_MOLECULARMATCH_CRITERIA_UNMET_SIZES);
                LOGGER.warn(keysCriteriaUnmet);
            }

            criteriaUnmetList.add(ImmutableMolecularMatchCriteriaUnmet.builder()
                    .priority(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("priority").getAsString())
                    .compositeKey(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("compositeKey").getAsString())
                    .isNew(!criteriaUnmet.getAsJsonObject().has("isNew")
                            ? null
                            : criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("isNew").getAsString())
                    .generatedBy(!criteriaUnmet.getAsJsonObject().has("generatedBy")
                            ? null
                            : criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("generatedBy").getAsString())
                    .manualSuppress(!criteriaUnmet.getAsJsonObject().has("manualSuppress")
                            ? null
                            : criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("manualSuppress").getAsString())
                    .generatedByTerm(!criteriaUnmet.getAsJsonObject().has("generatedByTerm")
                            ? null
                            : criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("generatedByTerm").getAsString())
                    .suppress(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("suppress").getAsString())
                    .filterType(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("filterType").getAsString())
                    .term(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("term").getAsString())
                    .primary(!criteriaUnmet.getAsJsonObject().has("primary")
                            ? null
                            : criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("primary").getAsString())
                    .facet(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("facet").getAsString())
                    .valid(!criteriaUnmet.getAsJsonObject().has("valid")
                            ? null
                            : criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("valid").getAsString())
                    .custom(!criteriaUnmet.getAsJsonObject().has("custom")
                            ? null
                            : criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("custom").getAsString())
                    .transcript(!criteriaUnmet.getAsJsonObject().has("transcript")
                            ? null
                            : criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("transcript").getAsString())
                    .build());
        }
        return criteriaUnmetList;
    }

    @NotNull
    private static List<MolecularMatchPrevalence> createPrevalence(@NotNull JsonArray arrayPrevelance) {
        List<MolecularMatchPrevalence> prevalenceList = Lists.newArrayList();

        for (JsonElement prevalence : arrayPrevelance) {
            Set<String> keysPrevalence = prevalence.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_PREVALENCE_SIZES.contains(keysPrevalence.size())) {
                LOGGER.warn("Found {} in molecular match prevalence rather than the expected {}",
                        keysPrevalence.size(),
                        EXPECTED_MOLECULARMATCH_PREVALENCE_SIZES);
                LOGGER.warn(keysPrevalence);
            }

            prevalenceList.add(ImmutableMolecularMatchPrevalence.builder()
                    .count(prevalence.getAsJsonObject().getAsJsonPrimitive("count").getAsString())
                    .percent(prevalence.getAsJsonObject().getAsJsonPrimitive("percent").getAsString())
                    .studyId(prevalence.getAsJsonObject().getAsJsonPrimitive("studyId").getAsString())
                    .samples(prevalence.getAsJsonObject().getAsJsonPrimitive("samples").getAsString())
                    .molecular(!prevalence.getAsJsonObject().has("molecular")
                            ? null
                            : prevalence.getAsJsonObject().getAsJsonPrimitive("molecular").getAsString())
                    .condition(!prevalence.getAsJsonObject().has("condition")
                            ? null
                            : prevalence.getAsJsonObject().getAsJsonPrimitive("condition").getAsString())
                    .build());
        }
        return prevalenceList;
    }

    @NotNull
    private static List<MolecularMatchMutations> createMutations(@NotNull JsonArray arrayMutations) {
        List<MolecularMatchMutations> mutationList = Lists.newArrayList();

        for (JsonElement mutation : arrayMutations) {
            Set<String> keysMutations = mutation.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_MUTATIONS_SIZES.contains(keysMutations.size())) {
                LOGGER.warn("Found {} in molecular match mutations rather than the expected {}",
                        keysMutations.size(),
                        EXPECTED_MOLECULARMATCH_MUTATIONS_SIZES);
                LOGGER.warn(keysMutations);
            }

            mutationList.add(ImmutableMolecularMatchMutations.builder()
                    .transcriptConsequence(mutation.getAsJsonObject().get("transcriptConsequence") == null
                            ? null
                            : createTranscriptConsequence(mutation.getAsJsonObject().getAsJsonArray("transcriptConsequence")))
                    .longestTranscript(!mutation.getAsJsonObject().has("longestTranscript")
                            ? null
                            : mutation.getAsJsonObject().getAsJsonPrimitive("longestTranscript").getAsString())
                    .parents(createParents(mutation.getAsJsonObject().getAsJsonArray("parents")))
                    .wgsaData(!mutation.getAsJsonObject().has("wgsaData")
                            ? null
                            : createWgsaData(mutation.getAsJsonObject().getAsJsonObject("wgsaData")))
                    .wgsaMap(!mutation.getAsJsonObject().has("wgsaMap")
                            ? null
                            : createWgsaMap(mutation.getAsJsonObject().getAsJsonArray("wgsaMap")))
                    .exonsInfo(!mutation.getAsJsonObject().has("exonsInfo")
                            ? null
                            : createExonsInfo(mutation.getAsJsonObject().getAsJsonObject("exonsInfo")))
                    .fusionData(!mutation.getAsJsonObject().has("fusionData")
                            ? null
                            : createFusionData(mutation.getAsJsonObject().getAsJsonArray("fusionData")))
                    .transcriptRecognized(!mutation.getAsJsonObject().has("transcriptRecognized")
                            ? null
                            : mutation.getAsJsonObject().get("transcriptRecognized").getAsString())
                    .description(mutation.getAsJsonObject().getAsJsonPrimitive("description").getAsString())
                    .mutationType(toStringList(mutation.getAsJsonObject().getAsJsonArray("mutation_type")))
                    .src(mutation.getAsJsonObject().getAsJsonPrimitive("_src").getAsString())
                    .sources(toStringList(mutation.getAsJsonObject().getAsJsonArray("sources")))
                    .synonyms(toStringList(mutation.getAsJsonObject().getAsJsonArray("synonyms")))
                    .gRch37Location(createGRCH37Location(mutation.getAsJsonObject().getAsJsonArray("GRCh37_location")))
                    .uniprotTranscript(!mutation.getAsJsonObject().has("uniprotTranscript")
                            ? null
                            : mutation.getAsJsonObject().getAsJsonPrimitive("uniprotTranscript").getAsString())
                    .geneSymbol(mutation.getAsJsonObject().getAsJsonPrimitive("geneSymbol").getAsString())
                    .pathology(toStringList(mutation.getAsJsonObject().getAsJsonArray("pathology")))
                    .transcript(!mutation.getAsJsonObject().has("transcript")
                            ? null
                            : mutation.getAsJsonObject().getAsJsonPrimitive("transcript").getAsString())
                    .id(mutation.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .cDNA(toStringList(mutation.getAsJsonObject().getAsJsonArray("cdna")))
                    .name(mutation.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .build());
        }
        return mutationList;
    }

    @NotNull
    private static List<MolecularMatchFusionData> createFusionData(@NotNull JsonArray arrayFusionData) {
        List<MolecularMatchFusionData> fusionDataList = Lists.newArrayList();
        for (JsonElement fusionData : arrayFusionData.getAsJsonArray()) {
            Set<String> keysFusionData = fusionData.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_FUSIONDATA_SIZES.contains(keysFusionData.size())) {
                LOGGER.warn("Found {} in molecular match fusion data rather than the expected {}",
                        keysFusionData.size(),
                        EXPECTED_MOLECULARMATCH_FUSIONDATA_SIZES);
                LOGGER.warn(keysFusionData);
            }

            fusionDataList.add(ImmutableMolecularMatchFusionData.builder()
                    .Bgreg(!fusionData.getAsJsonObject().has("Bgreg")
                            ? null
                            : createBreg(fusionData.getAsJsonObject().getAsJsonArray("Bgreg")))
                    .Bchr(!fusionData.getAsJsonObject().has("Bchr")
                            ? null
                            : toStringList(fusionData.getAsJsonObject().getAsJsonArray("Bchr")))
                    .synonym(!fusionData.getAsJsonObject().has("synonym")
                            ? null
                            : fusionData.getAsJsonObject().getAsJsonPrimitive("synonym").getAsString())
                    .Agene(!fusionData.getAsJsonObject().has("Agene")
                            ? null
                            : toStringList(fusionData.getAsJsonObject().getAsJsonArray("Agene")))
                    .Btx(!fusionData.getAsJsonObject().has("Btx")
                            ? null
                            : toStringList(fusionData.getAsJsonObject().getAsJsonArray("Btx")))
                    .Achr(!fusionData.getAsJsonObject().has("Achr")
                            ? null
                            : toStringList(fusionData.getAsJsonObject().getAsJsonArray("Achr")))
                    .ins(!fusionData.getAsJsonObject().has("ins")
                            ? null
                            : toStringList(fusionData.getAsJsonObject().getAsJsonArray("ins")))
                    .source(!fusionData.getAsJsonObject().has("source")
                            ? null
                            : fusionData.getAsJsonObject().getAsJsonPrimitive("source").getAsString())
                    .Agreg(!fusionData.getAsJsonObject().has("Agreg")
                            ? null
                            : createAreg(fusionData.getAsJsonObject().getAsJsonArray("Agreg")))
                    .Bgene(!fusionData.getAsJsonObject().has("Bgene")
                            ? null
                            : toStringList(fusionData.getAsJsonObject().getAsJsonArray("Bgene")))
                    .Acoord(!fusionData.getAsJsonObject().has("Acoord")
                            ? null
                            : toStringList(fusionData.getAsJsonObject().getAsJsonArray("Acoord")))
                    .Bori(!fusionData.getAsJsonObject().has("Bori")
                            ? null
                            : toStringList(fusionData.getAsJsonObject().getAsJsonArray("Bori")))
                    .Aband(!fusionData.getAsJsonObject().has("Aband")
                            ? null
                            : toStringList(fusionData.getAsJsonObject().getAsJsonArray("Aband")))
                    .Bband(!fusionData.getAsJsonObject().has("Bband")
                            ? null
                            : toStringList(fusionData.getAsJsonObject().getAsJsonArray("Bband")))
                    .Aori(!fusionData.getAsJsonObject().has("Aori")
                            ? null
                            : toStringList(fusionData.getAsJsonObject().getAsJsonArray("Aori")))
                    .Atx(!fusionData.getAsJsonObject().has("Atx")
                            ? null
                            : toStringList(fusionData.getAsJsonObject().getAsJsonArray("Atx")))
                    .Bcoord(!fusionData.getAsJsonObject().has("Bcoord")
                            ? null
                            : toStringList(fusionData.getAsJsonObject().getAsJsonArray("Bcoord")))
                    .Paper(!fusionData.getAsJsonObject().has("Paper")
                            ? null
                            : fusionData.getAsJsonObject().getAsJsonPrimitive("Paper").getAsString())
                    .build());
        }
        return fusionDataList;
    }

    @NotNull
    private static List<MolecularMatchAreg> createAreg(@NotNull JsonArray arrayAreg) {
        List<MolecularMatchAreg> aregList = Lists.newArrayList();
        for (JsonElement areg : arrayAreg.getAsJsonArray()) {
            Set<String> keysAreg = areg.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_AREG_SIZES.contains(keysAreg.size())) {
                LOGGER.warn("Found {} in molecular match areg rather than the expected {}",
                        keysAreg.size(),
                        EXPECTED_MOLECULARMATCH_AREG_SIZES);
                LOGGER.warn(keysAreg);
            }

            aregList.add(ImmutableMolecularMatchAreg.builder()
                    .type(areg.getAsJsonObject().getAsJsonPrimitive("type").getAsString())
                    .num(areg.getAsJsonObject().getAsJsonPrimitive("type").getAsString())
                    .build());
        }
        return aregList;
    }

    @NotNull
    private static List<MolecularMatchBreg> createBreg(@NotNull JsonArray arrayBreg) {
        List<MolecularMatchBreg> bregList = Lists.newArrayList();
        for (JsonElement breg : arrayBreg.getAsJsonArray()) {
            Set<String> keysBreg = breg.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_BREG_SIZES.contains(keysBreg.size())) {
                LOGGER.warn("Found {} in molecular match breg rather than the expected {}",
                        keysBreg.size(),
                        EXPECTED_MOLECULARMATCH_BREG_SIZES);
                LOGGER.warn(keysBreg);
            }

            bregList.add(ImmutableMolecularMatchBreg.builder()
                    .type(breg.getAsJsonObject().getAsJsonPrimitive("type").getAsString())
                    .num(breg.getAsJsonObject().getAsJsonPrimitive("num").getAsString())
                    .build());
        }
        return bregList;
    }

    @NotNull
    private static MolecularMatchExonInfo createExonsInfo(@NotNull JsonObject objectExonsInfo) {
        Set<String> keysExonsInfo = objectExonsInfo.keySet();
        if (!EXPECTED_MOLECULARMATCH_EXONSINFO_SIZES.contains(keysExonsInfo.size())) {
            LOGGER.warn("Found {} in molecular match exon info rather than the expected {}",
                    keysExonsInfo.size(),
                    EXPECTED_MOLECULARMATCH_EXONSINFO_SIZES);
            LOGGER.warn(keysExonsInfo);
        }

        return ImmutableMolecularMatchExonInfo.builder()
                .exonBoundaries(createExonBoundaries(objectExonsInfo.getAsJsonObject("exonBoundaries")))
                .txStart(!objectExonsInfo.has("txStart") ? null : objectExonsInfo.getAsJsonPrimitive("txStart").getAsString())
                .cdsEnd(!objectExonsInfo.has("cdsEnd") ? null : objectExonsInfo.getAsJsonPrimitive("cdsEnd").getAsString())
                .chr(objectExonsInfo.getAsJsonPrimitive("chr").getAsString())
                .cdsStart(!objectExonsInfo.has("cdsEnd") ? null : objectExonsInfo.getAsJsonPrimitive("cdsEnd").getAsString())
                .transcript(objectExonsInfo.getAsJsonPrimitive("transcript").getAsString())
                .txEnd(!objectExonsInfo.has("txEnd") ? null : objectExonsInfo.getAsJsonPrimitive("txEnd").getAsString())
                .build();
    }

    @NotNull
    private static MolecularMatchPositions createMolecularPositions(@NotNull JsonObject objectPositions) {
        Set<String> keysPositions = objectPositions.keySet();
        if (!EXPECTED_MOLECULARMATCH_POSITIONS_SIZES.contains(keysPositions.size())) {
            LOGGER.warn("Found {} in molecular match positions rather than the expected {}",
                    keysPositions.size(),
                    EXPECTED_MOLECULARMATCH_POSITIONS_SIZES);
            LOGGER.warn(keysPositions);
        }

        return ImmutableMolecularMatchPositions.builder()
                .start(objectPositions.getAsJsonPrimitive("start").getAsString())
                .stop(objectPositions.getAsJsonPrimitive("stop").getAsString())
                .build();
    }

    @NotNull
    private static MolecularMatchExonsBoundries createExonBoundaries(@NotNull JsonObject objectExonsBoundaries) {
        Set<String> keysExonsBoundaries = objectExonsBoundaries.keySet();
        if (!EXPECTED_MOLECULARMATCH_EXONSBOUNDARIES_SIZES.contains(keysExonsBoundaries.size())) {
            LOGGER.warn("Found {} in molecular match exons boundaries rather than the expected {}",
                    keysExonsBoundaries.size(),
                    EXPECTED_MOLECULARMATCH_EXONSBOUNDARIES_SIZES);
            LOGGER.warn(keysExonsBoundaries);
        }

        return ImmutableMolecularMatchExonsBoundries.builder()
                .exon1(!objectExonsBoundaries.has("1") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("1")))
                .exon2(!objectExonsBoundaries.has("2") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("2")))
                .exon3(!objectExonsBoundaries.has("3") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("3")))
                .exon4(!objectExonsBoundaries.has("4") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("4")))
                .exon5(!objectExonsBoundaries.has("5") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("5")))
                .exon6(!objectExonsBoundaries.has("6") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("6")))
                .exon7(!objectExonsBoundaries.has("7") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("7")))
                .exon8(!objectExonsBoundaries.has("8") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("8")))
                .exon9(!objectExonsBoundaries.has("9") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("9")))
                .exon10(!objectExonsBoundaries.has("10") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("10")))
                .exon11(!objectExonsBoundaries.has("11") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("11")))
                .exon12(!objectExonsBoundaries.has("12") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("12")))
                .exon13(!objectExonsBoundaries.has("13") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("13")))
                .exon14(!objectExonsBoundaries.has("14") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("14")))
                .exon15(!objectExonsBoundaries.has("15") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("15")))
                .exon16(!objectExonsBoundaries.has("16") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("16")))
                .exon17(!objectExonsBoundaries.has("17") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("17")))
                .exon18(!objectExonsBoundaries.has("18") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("18")))
                .exon19(!objectExonsBoundaries.has("19") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("19")))
                .exon20(!objectExonsBoundaries.has("20") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("20")))
                .exon21(!objectExonsBoundaries.has("21") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("21")))
                .exon22(!objectExonsBoundaries.has("22") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("22")))
                .exon23(!objectExonsBoundaries.has("23") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("23")))
                .exon24(!objectExonsBoundaries.has("24") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("24")))
                .exon25(!objectExonsBoundaries.has("25") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("25")))
                .exon26(!objectExonsBoundaries.has("26") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("26")))
                .exon27(!objectExonsBoundaries.has("27") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("27")))
                .exon28(!objectExonsBoundaries.has("28") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("28")))
                .exon29(!objectExonsBoundaries.has("29") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("29")))
                .exon30(!objectExonsBoundaries.has("30") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("30")))
                .exon31(!objectExonsBoundaries.has("31") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("31")))
                .exon32(!objectExonsBoundaries.has("32") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("32")))
                .exon33(!objectExonsBoundaries.has("33") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("33")))
                .exon34(!objectExonsBoundaries.has("34") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("34")))
                .exon35(!objectExonsBoundaries.has("35") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("35")))
                .exon36(!objectExonsBoundaries.has("36") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("36")))
                .exon37(!objectExonsBoundaries.has("37") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("37")))
                .exon38(!objectExonsBoundaries.has("38") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("38")))
                .exon39(!objectExonsBoundaries.has("39") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("39")))
                .exon40(!objectExonsBoundaries.has("40") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("40")))
                .exon41(!objectExonsBoundaries.has("41") ? null : createMolecularPositions(objectExonsBoundaries.getAsJsonObject("41")))
                .build();
    }

    @NotNull
    private static List<MolecularMatchWGSaMap> createWgsaMap(@NotNull JsonArray objectWgsaMap) {
        List<MolecularMatchWGSaMap> molecularMatchWGSaMapList = Lists.newArrayList();
        for (JsonElement wgsDataMap : objectWgsaMap.getAsJsonArray()) {
            Set<String> keysWgsaMap = wgsDataMap.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_WGSAMAP_SIZES.contains(keysWgsaMap.size())) {
                LOGGER.warn("Found {} in molecular match wgsa map rather than the expected {}",
                        keysWgsaMap.size(),
                        EXPECTED_MOLECULARMATCH_WGSAMAP_SIZES);
                LOGGER.warn(keysWgsaMap);
            }

            molecularMatchWGSaMapList.add(ImmutableMolecularMatchWGSaMap.builder()
                    .AA(!wgsDataMap.getAsJsonObject().has("AA")
                            ? null
                            : wgsDataMap.getAsJsonObject().getAsJsonPrimitive("AA").getAsString())
                    .name(wgsDataMap.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .GRCh37_Chr_Start_Ref_Alt(wgsDataMap.getAsJsonObject().getAsJsonPrimitive("GRCh37_Chr_Start_Ref_Alt").getAsString())
                    .Synonyms(toStringList(wgsDataMap.getAsJsonObject().getAsJsonArray("Synonyms")))
                    .ProtCoords(toStringList(wgsDataMap.getAsJsonObject().getAsJsonArray("ProtCoords")))
                    .NucleotideChange(wgsDataMap.getAsJsonObject().getAsJsonPrimitive("NucleotideChange").getAsString())
                    .Exon(!wgsDataMap.getAsJsonObject().has("Exon")
                            ? null
                            : wgsDataMap.getAsJsonObject().getAsJsonPrimitive("Exon").getAsString())
                    .Gene(wgsDataMap.getAsJsonObject().getAsJsonPrimitive("Gene").getAsString())
                    .Transcript(wgsDataMap.getAsJsonObject().getAsJsonPrimitive("Transcript").getAsString())
                    .build());
        }
        return molecularMatchWGSaMapList;
    }

    @NotNull
    private static List<MolecularMatchWGSadataLocation> createWgsaData(@NotNull JsonObject objectWgsaData) {
        List<MolecularMatchWGSadataLocation> molecularMatchWGSadataLocationList = Lists.newArrayList();
        Set<String> keysWgsaData = objectWgsaData.keySet();
        if (!EXPECTED_MOLECULARMATCH_WGSADATA_SIZES.contains(keysWgsaData.size())) {
            LOGGER.warn("Found {} in molecular match wgsa data rather than the expected {}",
                    keysWgsaData.size(),
                    EXPECTED_MOLECULARMATCH_WGSADATA_SIZES);
            LOGGER.warn(keysWgsaData);
        }

        for (JsonElement wgsDataLocation : objectWgsaData.get("locations").getAsJsonArray()) {
            Set<String> keysWgsaDataLocation = wgsDataLocation.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_WGSADATA_LOCATION_SIZES.contains(keysWgsaDataLocation.size())) {
                LOGGER.warn("Found {} in molecular match wgsa data locations rather than the expected {}",
                        keysWgsaDataLocation.size(),
                        EXPECTED_MOLECULARMATCH_WGSADATA_LOCATION_SIZES);
                LOGGER.warn(keysWgsaDataLocation);
            }

            molecularMatchWGSadataLocationList.add(ImmutableMolecularMatchWGSadataLocation.builder()
                    .ExonicFunc(!wgsDataLocation.getAsJsonObject().has("ExonicFunc")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("ExonicFunc").getAsString())
                    .dbSNP(!wgsDataLocation.getAsJsonObject().has("dbSNP")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("dbSNP").getAsString())
                    .ClinVar_DIS(!wgsDataLocation.getAsJsonObject().has("ClinVar_DIS")
                            ? null
                            : toStringList(wgsDataLocation.getAsJsonObject().getAsJsonArray("ClinVar_DIS")))
                    .ClinVar_SIG(!wgsDataLocation.getAsJsonObject().has("ClinVar_SIG")
                            ? null
                            : toStringList(wgsDataLocation.getAsJsonObject().getAsJsonArray("ClinVar_SIG")))
                    .ClinVar_STATUS(!wgsDataLocation.getAsJsonObject().has("ClinVar_STATUS")
                            ? null
                            : toStringList(wgsDataLocation.getAsJsonObject().getAsJsonArray("ClinVar_STATUS")))
                    .ClinVar_DBID(!wgsDataLocation.getAsJsonObject().has("ClinVar_DBID")
                            ? null
                            : toStringList(wgsDataLocation.getAsJsonObject().getAsJsonArray("ClinVar_DBID")))
                    .ExAC_NFE(!wgsDataLocation.getAsJsonObject().has("ExAC_NFE")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("ExAC_NFE").getAsString())
                    .ExAC_FIN(!wgsDataLocation.getAsJsonObject().has("ExAC_FIN")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("ExAC_FIN").getAsString())
                    .G1000_ALL(!wgsDataLocation.getAsJsonObject().has("1000G_ALL")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("1000G_ALL").getAsString())
                    .G1000_SAS(!wgsDataLocation.getAsJsonObject().has("1000G_SAS")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("1000G_SAS").getAsString())
                    .G1000_EAS(!wgsDataLocation.getAsJsonObject().has("1000G_EAS")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("1000G_EAS").getAsString())
                    .G1000_AFR(!wgsDataLocation.getAsJsonObject().has("1000G_AFR")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("1000G_AFR").getAsString())
                    .ExAC_SAS(!wgsDataLocation.getAsJsonObject().has("ExAC_SAS")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("ExAC_SAS").getAsString())
                    .ExAC_EAS(!wgsDataLocation.getAsJsonObject().has("ExAC_EAS")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("ExAC_EAS").getAsString())
                    .ExAC_AMR(!wgsDataLocation.getAsJsonObject().has("ExAC_AMR")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("ExAC_AMR").getAsString())
                    .ExAC_AFR(!wgsDataLocation.getAsJsonObject().has("ExAC_AFR")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("ExAC_AFR").getAsString())
                    .ExAC_Freq(!wgsDataLocation.getAsJsonObject().has("ExAC_Freq")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("ExAC_Freq").getAsString())
                    .End(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("End").getAsString())
                    .Start(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("Start").getAsString())
                    .SiPhy_29way_logOdds(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("SiPhy_29way_logOdds").getAsString())
                    .FullAA(toStringList(wgsDataLocation.getAsJsonObject().getAsJsonArray("FullAA")))
                    .Ref(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("Ref").getAsString())
                    .GERP_RS(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("GERP++_RS").getAsString())
                    .FATHMM(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("FATHMM").getAsString())
                    .NucleotideChange(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("NucleotideChange").getAsString())
                    .phyloP100way_vertebrate(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("phyloP100way_vertebrate").getAsString())
                    .Func(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("phyloP100way_vertebrate").getAsString())
                    .GWAS_PUBMED(!wgsDataLocation.getAsJsonObject().has("GWAS_PUBMED")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("GWAS_PUBMED").getAsString())
                    .Transcript(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("Transcript").getAsString())
                    .ESP6500si_AA(!wgsDataLocation.getAsJsonObject().has("ESP6500si_AA")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("ESP6500si_AA").getAsString())
                    .ESP6500si_EA(!wgsDataLocation.getAsJsonObject().has("ESP6500si_EA")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("ESP6500si_EA").getAsString())
                    .G1000_EUR(!wgsDataLocation.getAsJsonObject().has("1000G_EUR")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("1000G_EUR").getAsString())
                    .G1000_AMR(!wgsDataLocation.getAsJsonObject().has("1000G_AMR")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("1000G_AMR").getAsString())
                    .Chr_Start_Ref_Alt(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("Chr_Start_Ref_Alt").getAsString())
                    .AA(!wgsDataLocation.getAsJsonObject().has("AA")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("AA").getAsString())
                    .PopFreqMax(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("PopFreqMax").getAsString())
                    .FATHMM_Pred(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("FATHMM_Pred").getAsString())
                    .wgRna(!wgsDataLocation.getAsJsonObject().has("wgRna")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("wgRna").getAsString())
                    .Gene(toStringList(wgsDataLocation.getAsJsonObject().getAsJsonArray("Gene")))
                    .phyloP46way_placental(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("phyloP46way_placental").getAsString())
                    .key(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("_key").getAsString())
                    .targetScanS(!wgsDataLocation.getAsJsonObject().has("targetScanS")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("targetScanS").getAsString())
                    .Chr(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("Chr").getAsString())
                    .COSMIC_ID(!wgsDataLocation.getAsJsonObject().has("COSMIC_ID")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("COSMIC_ID").getAsString())
                    .alt(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("Alt").getAsString())
                    .GWAS_DIS(!wgsDataLocation.getAsJsonObject().has("GWAS_DIS")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("GWAS_DIS").getAsString())
                    .GWAS_SNP(!wgsDataLocation.getAsJsonObject().has("GWAS_SNP")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("GWAS_SNP").getAsString())
                    .build());
        }
        return molecularMatchWGSadataLocationList;
    }

    @NotNull
    private static List<MolecularMatchGRch37Location> createGRCH37Location(@NotNull JsonArray arrayLocation) {
        List<MolecularMatchGRch37Location> grch37LocationList = Lists.newArrayList();

        for (JsonElement location : arrayLocation) {
            Set<String> keysLocation = location.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_LOCATIONGRCH37_SIZES.contains(keysLocation.size())) {
                LOGGER.warn("Found {} in molecular match grch37 location rather than the expected {}",
                        keysLocation.size(),
                        EXPECTED_MOLECULARMATCH_LOCATIONGRCH37_SIZES);
                LOGGER.warn(keysLocation);
            }

            grch37LocationList.add(ImmutableMolecularMatchGRch37Location.builder()
                    .compositeKey(location.getAsJsonObject().getAsJsonPrimitive("compositeKey").getAsString())
                    .ref(location.getAsJsonObject().get("ref").isJsonNull()
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("ref").getAsString())
                    .stop(location.getAsJsonObject().get("stop").isJsonNull()
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("stop").getAsString())
                    .start(location.getAsJsonObject().get("start").isJsonNull()
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("start").getAsString())
                    .chr(location.getAsJsonObject().get("chr").isJsonNull()
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("chr").getAsString())
                    .alt(location.getAsJsonObject().get("alt").isJsonNull()
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("alt").getAsString())
                    .validated(location.getAsJsonObject().getAsJsonPrimitive("validated").getAsString())
                    .transcriptConsequences(createConsequencesGRCH37(location.getAsJsonObject().getAsJsonArray("transcript_consequences")))
                    .strand(location.getAsJsonObject().getAsJsonPrimitive("strand").getAsString())
                    .build());
        }
        return grch37LocationList;
    }

    @NotNull
    private static List<MolecularMatchTranscriptConsequencesGRCH37> createConsequencesGRCH37(
            @NotNull JsonArray arrayTranscriptConsequence) {
        List<MolecularMatchTranscriptConsequencesGRCH37> transcriptConsequencesGRCH37List = Lists.newArrayList();
        for (JsonElement transcriptConsequences : arrayTranscriptConsequence) {
            Set<String> keysTranscriptConsequences = transcriptConsequences.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_TRANSCRIPT_CONSEQUENCES__GRCH37_SIZES.contains(keysTranscriptConsequences.size())) {
                LOGGER.warn("Found {} in molecular match transcript consequences grch37 rather than the expected {}",
                        keysTranscriptConsequences.size(),
                        EXPECTED_MOLECULARMATCH_TRANSCRIPT_CONSEQUENCES__GRCH37_SIZES);
                LOGGER.warn(keysTranscriptConsequences);
            }

            transcriptConsequencesGRCH37List.add(ImmutableMolecularMatchTranscriptConsequencesGRCH37.builder()
                    .aminoAcidChange(transcriptConsequences.getAsJsonObject().get("amino_acid_change").isJsonNull()
                            ? null
                            : transcriptConsequences.getAsJsonObject().getAsJsonPrimitive("amino_acid_change").getAsString())
                    .txSites(toStringList(transcriptConsequences.getAsJsonObject().getAsJsonArray("txSites")))
                    .exonNumber(transcriptConsequences.getAsJsonObject().get("exonNumber").isJsonNull()
                            ? null
                            : createArrayExonNumber(transcriptConsequences))
                    .intronNumber(transcriptConsequences.getAsJsonObject().get("intronNumber").isJsonNull()
                            ? null
                            : transcriptConsequences.getAsJsonObject().getAsJsonPrimitive("intronNumber").getAsString())
                    .transcript(transcriptConsequences.getAsJsonObject().getAsJsonPrimitive("transcript").getAsString())
                    .cdna(transcriptConsequences.getAsJsonObject().get("cdna").isJsonNull()
                            ? null
                            : transcriptConsequences.getAsJsonObject().getAsJsonPrimitive("cdna").getAsString())
                    .build());
        }
        return transcriptConsequencesGRCH37List;
    }

    @NotNull
    private static List<MolecularMatchTranscriptConsequence> createTranscriptConsequence(@NotNull JsonArray arrayTranscriptConsequence) {
        List<MolecularMatchTranscriptConsequence> transcriptConsequenceList = Lists.newArrayList();

        for (JsonElement transcriptConsequence : arrayTranscriptConsequence) {
            Set<String> keysTranscriptConsequence = transcriptConsequence.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_TRANSCRIPT_CONSEQUENCES_SIZES.contains(keysTranscriptConsequence.size())) {
                LOGGER.warn("Found {} in molecular match transcript consequence rather than the expected {}",
                        keysTranscriptConsequence.size(),
                        EXPECTED_MOLECULARMATCH_TRANSCRIPT_CONSEQUENCES_SIZES);
                LOGGER.warn(keysTranscriptConsequence);
            }

            transcriptConsequenceList.add(ImmutableMolecularMatchTranscriptConsequence.builder()
                    .aminoAcidChange(
                            !transcriptConsequence.getAsJsonObject().has("amino_acid_change") || transcriptConsequence.getAsJsonObject()
                                    .get("amino_acid_change")
                                    .isJsonNull() ? null : transcriptConsequence.getAsJsonObject().get("amino_acid_change").getAsString())
                    .compositeKey(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("compositeKey").getAsString())
                    .intronNumber(transcriptConsequence.getAsJsonObject().get("intronNumber").isJsonNull()
                            ? null
                            : transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("intronNumber").getAsString())
                    .exonNumber(transcriptConsequence.getAsJsonObject().get("exonNumber").isJsonNull()
                            ? null
                            : createArrayExonNumber(transcriptConsequence))
                    .suppress(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("suppress").getAsString())
                    .stop(!transcriptConsequence.getAsJsonObject().has("stop")
                            ? null
                            : transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("stop").getAsString())
                    .custom(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("custom").getAsString())
                    .start(!transcriptConsequence.getAsJsonObject().has("start")
                            ? null
                            : transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("start").getAsString())
                    .chr(!transcriptConsequence.getAsJsonObject().has("chr")
                            ? null
                            : transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("chr").getAsString())
                    .strand(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("strand").getAsString())
                    .validated(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("validated").getAsString())
                    .transcript(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("transcript").getAsString())
                    .cdna(!transcriptConsequence.getAsJsonObject().has("cdna")
                            ? null
                            : transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("cdna").getAsString())
                    .referenceGenome(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("referenceGenome").getAsString())
                    .ref(!transcriptConsequence.getAsJsonObject().has("ref")
                            ? null
                            : transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("ref").getAsString())
                    .alt(!transcriptConsequence.getAsJsonObject().has("alt")
                            ? null
                            : transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("alt").getAsString())
                    .build());
        }
        return transcriptConsequenceList;
    }
}
