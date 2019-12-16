package com.hartwig.hmftools.vicc.reader;

import static com.hartwig.hmftools.vicc.reader.JsonFunctions.nullableString;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.optionalJsonArray;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.optionalJsonObject;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.optionalNullableString;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.string;
import static com.hartwig.hmftools.vicc.reader.JsonFunctions.stringList;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.vicc.datamodel.civic.Civic;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicAvatars;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicClinicalTrial;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicCoordinates;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicDisease;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicDrug;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicEvidenceItem;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicLastCommentedOn;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicLastModified;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicLastReviewed;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicLifecycleActions;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicOrganization;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicProfileImage;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicProvisionalValue;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicPublicationDate;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicSource;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicUser;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicVariant;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicVariantGroup;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicVariantType;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivic;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicAvatars;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicClinicalTrial;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicCoordinates;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicDisease;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicDrug;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicEvidenceItem;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicLastCommentedOn;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicLastModified;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicLastReviewed;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicLifecycleActions;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicOrganization;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicProfileImage;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicProvisionalValue;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicPublicationDate;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicSource;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicUser;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicVariant;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicVariantGroup;
import com.hartwig.hmftools.vicc.datamodel.civic.ImmutableCivicVariantType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class CivicObjectFactory {

    private static final Logger LOGGER = LogManager.getLogger(CivicObjectFactory.class);

    private static final List<Integer> EXPECTED_CIVIC_AVATARS_SIZES = Lists.newArrayList(4);
    private static final List<Integer> EXPECTED_CIVIC_DISEASES_SIZES = Lists.newArrayList(5);
    private static final List<Integer> EXPECTED_CIVIC_DRUGS_SIZES = Lists.newArrayList(3);
    private static final List<Integer> EXPECTED_CIVIC_EVIDENCE_ITEMS_SIZES = Lists.newArrayList(17);
    private static final List<Integer> EXPECTED_CIVIC_LAST_COMMENTED_ON_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_CIVIC_LAST_MODIFIED_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_CIVIC_LAST_REVIEWED_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_CIVIC_LIFECYCLE_ACTIONS_SIZES = Lists.newArrayList(0, 3);
    private static final List<Integer> EXPECTED_CIVIC_ORGANIZATION_SIZES = Lists.newArrayList(0, 5);
    private static final List<Integer> EXPECTED_CIVIC_PROFILE_IMAGE_SIZES = Lists.newArrayList(5);
    private static final List<Integer> EXPECTED_CIVIC_PUBLICATIONS_DATE_SIZES = Lists.newArrayList(0, 1, 2, 3);
    private static final List<Integer> EXPECTED_CIVIC_SOURCE_SIZES = Lists.newArrayList(13);
    private static final List<Integer> EXPECTED_CIVIC_USER_SIZES = Lists.newArrayList(21);
    private static final List<Integer> EXPECTED_CIVIC_VARIANTTYPES_SIZES = Lists.newArrayList(6);
    private static final List<Integer> EXPECTED_CIVIC_CLINICALTRIAL_SIZES = Lists.newArrayList(4);

    private CivicObjectFactory() {
    }

    @NotNull
    static Civic create(@NotNull JsonObject civicObject) {
        ViccDatamodelCheckerFactory.civicEntryChecker().check(civicObject);

        CivicCoordinates coordinates = createCoordinates(civicObject.getAsJsonObject("coordinates"));
        assert coordinates != null;

        return ImmutableCivic.builder()
                .entrezId(string(civicObject, "entrez_id"))
                .entrezName(string(civicObject, "entrez_name"))
                .name(string(civicObject, "name"))
                .type(string(civicObject, "type"))
                .coordinates(coordinates)
                .sources(createSources(civicObject.getAsJsonArray("sources")))
                .variantAliases(stringList(civicObject, "variant_aliases"))
                .variantGroups(createVariantGroups(optionalJsonArray(civicObject, "variant_groups")))
                .variantTypes(createVariantTypes(civicObject.getAsJsonArray("variant_types")))
                .hgvsExpressions(stringList(civicObject, "hgvs_expressions"))
                .evidenceItems(createEvidenceItems(civicObject.getAsJsonArray("evidence_items")))
                .assertions(stringList(civicObject, "assertions"))
                .civicActionabilityScore(nullableString(civicObject, "civic_actionability_score"))
                .clinvarEntries(stringList(civicObject, "clinvar_entries"))
                .alleleRegistryId(nullableString(civicObject, "allele_registry_id"))
                .provisionalValue(createProvisionalValue(civicObject.getAsJsonObject("provisional_values")))
                .lifecycleActions(createLifeCycleActions(civicObject.getAsJsonObject("lifecycle_actions")))
                .id(string(civicObject, "id"))
                .geneId(string(civicObject, "gene_id"))
                .description(string(civicObject, "description"))
                .build();
    }

    @Nullable
    private static CivicCoordinates createCoordinates(@Nullable JsonObject coordinatesObject) {
        if (coordinatesObject == null) {
            return null;

        }
        ViccDatamodelCheckerFactory.civicCoordinatesChecker().check(coordinatesObject);

        return ImmutableCivicCoordinates.builder()
                .chromosome(nullableString(coordinatesObject, "chromosome"))
                .start(nullableString(coordinatesObject, "start"))
                .stop(nullableString(coordinatesObject, "stop"))
                .referenceBases(nullableString(coordinatesObject, "reference_bases"))
                .variantBases(nullableString(coordinatesObject, "variant_bases"))
                .representativeTranscript(nullableString(coordinatesObject, "representative_transcript"))
                .ensemblVersion(nullableString(coordinatesObject, "ensembl_version"))
                .chromosome2(nullableString(coordinatesObject, "chromosome2"))
                .start2(nullableString(coordinatesObject, "start2"))
                .stop2(nullableString(coordinatesObject, "stop2"))
                .representativeTranscript2(nullableString(coordinatesObject, "representative_transcript2"))
                .build();
    }

    @NotNull
    private static List<CivicSource> createSources(@NotNull JsonArray sourceArray) {
        List<CivicSource> sourceList = Lists.newArrayList();
        ViccDatamodelChecker sourceChecker = ViccDatamodelCheckerFactory.civicSourceChecker();

        for (JsonElement sourceElement : sourceArray) {
            JsonObject sourceObject = sourceElement.getAsJsonObject();
            sourceChecker.check(sourceObject);

            sourceList.add(ImmutableCivicSource.builder()
                    .name(nullableString(sourceObject, "name"))
                    .status(string(sourceObject, "status"))
                    .openAccess(nullableString(sourceObject, "open_access"))
                    .journal(nullableString(sourceObject, "journal"))
                    .fullJournalTitle(nullableString(sourceObject, "full_journal_title"))
                    .citation(string(sourceObject, "citation"))
                    .pmcId(nullableString(sourceObject, "pmc_id"))
                    .sourceUrl(string(sourceObject, "source_url"))
                    .clinicalTrials(createCivicClinicalTrials(sourceObject.getAsJsonArray("clinical_trials")))
                    .pubmedId(string(sourceObject, "pubmed_id"))
                    .isReview(string(sourceObject, "is_review"))
                    .publicationDate(createPublicationDate(sourceObject.getAsJsonObject("publication_date")))
                    .id(string(sourceObject, "id"))
                    .build());
        }

        return sourceList;
    }

    @Nullable
    private static CivicProvisionalValue createProvisionalValue(@NotNull JsonObject provisionalValueObject) {
        ViccDatamodelCheckerFactory.civicProvisionalValueChecker().check(provisionalValueObject);

        JsonObject descriptionObject = optionalJsonObject(provisionalValueObject, "description");
        if (descriptionObject == null) {
            return null;
        }

        ViccDatamodelCheckerFactory.civicProvisionalValueDescriptionChecker().check(descriptionObject);

        return ImmutableCivicProvisionalValue.builder()
                .revisionId(optionalNullableString(descriptionObject, "revision_id"))
                .value(optionalNullableString(descriptionObject, "value"))
                .build();
    }

    @Nullable
    private static List<CivicVariantGroup> createVariantGroups(@Nullable JsonArray variantGroupArray) {
        if (variantGroupArray == null) {
            return null;
        }

        List<CivicVariantGroup> variantGroupList = Lists.newArrayList();
        ViccDatamodelChecker variantGroupChecker = ViccDatamodelCheckerFactory.civicVariantGroupChecker();
        for (JsonElement variantGroupElement : variantGroupArray) {
            JsonObject variantGroupObject = variantGroupElement.getAsJsonObject();
            variantGroupChecker.check(variantGroupObject);

            variantGroupList.add(ImmutableCivicVariantGroup.builder()
                    .name(string(variantGroupObject, "name"))
                    .type(string(variantGroupObject, "type"))
                    .description(string(variantGroupObject, "description"))
                    .variants(createVariants(variantGroupObject.getAsJsonArray("variants")))
                    .id(string(variantGroupObject, "id"))
                    .build());
        }
        return variantGroupList;
    }

    @NotNull
    private static List<CivicVariant> createVariants(@NotNull JsonArray variantArray) {
        List<CivicVariant> variantList = Lists.newArrayList();
        ViccDatamodelChecker variantChecker = ViccDatamodelCheckerFactory.civicVariantChecker();

        for (JsonElement variantElement : variantArray) {
            JsonObject variantObject = variantElement.getAsJsonObject();
            variantChecker.check(variantObject);

            variantList.add(ImmutableCivicVariant.builder()
                    .entrezId(string(variantObject, "entrez_id"))
                    .entrezName(string(variantObject, "entrez_name"))
                    .name(string(variantObject, "name"))
                    .type(string(variantObject, "type"))
                    .variantTypes(createVariantTypes(variantObject.getAsJsonArray("variant_types")))
                    .civicActionabilityScore(optionalNullableString(variantObject, "civic_actionability_score"))
                    .coordinates(createCoordinates(optionalJsonObject(variantObject, "createCoordinates")))
                    .id(string(variantObject, "id"))
                    .geneId(string(variantObject, "gene_id"))
                    .description(string(variantObject, "description"))
                    .build());
        }

        return variantList;
    }

    @NotNull
    private static List<CivicEvidenceItem> createEvidenceItems(@NotNull JsonArray evidenceItemsArray) {
        List<CivicEvidenceItem> evidenceItemsList = Lists.newArrayList();
        for (JsonElement evidenceItem : evidenceItemsArray) {
            Set<String> keysEvidenceItems = evidenceItem.getAsJsonObject().keySet();
            if (!EXPECTED_CIVIC_EVIDENCE_ITEMS_SIZES.contains(keysEvidenceItems.size())) {
                LOGGER.warn("Found {} in civic evidence items rather than the expected {}",
                        keysEvidenceItems.size(),
                        EXPECTED_CIVIC_EVIDENCE_ITEMS_SIZES);
                LOGGER.warn(keysEvidenceItems);
            }

            evidenceItemsList.add(ImmutableCivicEvidenceItem.builder()
                    .status(evidenceItem.getAsJsonObject().getAsJsonPrimitive("status").getAsString())
                    .rating(evidenceItem.getAsJsonObject().get("rating").isJsonNull()
                            ? null
                            : evidenceItem.getAsJsonObject().getAsJsonPrimitive("rating").getAsString())
                    .drugInteractionType(evidenceItem.getAsJsonObject().get("drug_interaction_type").isJsonNull()
                            ? null
                            : evidenceItem.getAsJsonObject().getAsJsonPrimitive("drug_interaction_type").getAsString())
                    .description(evidenceItem.getAsJsonObject().getAsJsonPrimitive("description").getAsString())
                    .openChangeCount(evidenceItem.getAsJsonObject().getAsJsonPrimitive("open_change_count").getAsString())
                    .evidenceType(evidenceItem.getAsJsonObject().getAsJsonPrimitive("evidence_type").getAsString())
                    .drugs(createDrugs(evidenceItem.getAsJsonObject().getAsJsonArray("drugs")))
                    .variantOrigin(evidenceItem.getAsJsonObject().get("variant_origin").isJsonNull()
                            ? null
                            : evidenceItem.getAsJsonObject().getAsJsonPrimitive("variant_origin").getAsString())
                    .disease(createDiseases(evidenceItem.getAsJsonObject().getAsJsonObject("disease")))
                    .source(createSource(evidenceItem.getAsJsonObject().getAsJsonObject("source")))
                    .evidenceDirection(evidenceItem.getAsJsonObject().get("evidence_direction").isJsonNull()
                            ? null
                            : evidenceItem.getAsJsonObject().getAsJsonPrimitive("evidence_direction").toString())
                    .variantId(!evidenceItem.getAsJsonObject().has("variant_id")
                            ? null
                            : evidenceItem.getAsJsonObject().getAsJsonPrimitive("variant_id").getAsString())
                    .clinicalSignificance(evidenceItem.getAsJsonObject().get("clinical_significance").isJsonNull()
                            ? null
                            : evidenceItem.getAsJsonObject().getAsJsonPrimitive("clinical_significance").getAsString())
                    .evidenceLevel(evidenceItem.getAsJsonObject().getAsJsonPrimitive("evidence_level").getAsString())
                    .type(evidenceItem.getAsJsonObject().getAsJsonPrimitive("type").getAsString())
                    .id(evidenceItem.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .name(evidenceItem.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .build());
        }
        return evidenceItemsList;
    }

    @NotNull
    private static CivicSource createSource(@NotNull JsonObject objectSource) {
        Set<String> keysSource = objectSource.keySet();
        if (!EXPECTED_CIVIC_SOURCE_SIZES.contains(keysSource.size())) {
            LOGGER.warn("Found {} in civic source rather than the expected {}", keysSource.size(), EXPECTED_CIVIC_SOURCE_SIZES);
            LOGGER.warn(keysSource);
        }

        return ImmutableCivicSource.builder()
                .status(objectSource.getAsJsonPrimitive("status").getAsString())
                .openAccess(objectSource.get("open_access").isJsonNull()
                        ? null
                        : objectSource.getAsJsonPrimitive("open_access").getAsString())
                .name(objectSource.get("name").isJsonNull() ? null : objectSource.getAsJsonPrimitive("name").getAsString())
                .journal(objectSource.get("journal").isJsonNull() ? null : objectSource.getAsJsonPrimitive("journal").getAsString())
                .citation(objectSource.getAsJsonPrimitive("citation").getAsString())
                .pmcId(objectSource.get("pmc_id").isJsonNull() ? null : objectSource.getAsJsonPrimitive("pmc_id").getAsString())
                .fullJournalTitle(objectSource.get("full_journal_title").isJsonNull()
                        ? null
                        : objectSource.getAsJsonPrimitive("full_journal_title").getAsString())
                .sourceUrl(objectSource.getAsJsonPrimitive("source_url").getAsString())
                .clinicalTrials(createCivicClinicalTrials(objectSource.getAsJsonArray("clinical_trials")))
                .pubmedId(objectSource.getAsJsonPrimitive("pubmed_id").getAsString())
                .isReview(objectSource.getAsJsonPrimitive("is_review").getAsString())
                .publicationDate(createPublicationDate(objectSource.getAsJsonObject("publication_date")))
                .id(objectSource.getAsJsonPrimitive("id").getAsString())
                .build();
    }

    @NotNull
    private static List<CivicClinicalTrial> createCivicClinicalTrials(@NotNull JsonArray arrayClinicalTrial) {
        List<CivicClinicalTrial> clinicalTrialList = Lists.newArrayList();
        for (JsonElement clinicalTrial : arrayClinicalTrial) {
            Set<String> keysClinicalTrials = clinicalTrial.getAsJsonObject().keySet();
            if (!EXPECTED_CIVIC_CLINICALTRIAL_SIZES.contains(keysClinicalTrials.size())) {
                LOGGER.warn("Found {} in civic clinical trials rather than the expected {}",
                        keysClinicalTrials.size(),
                        EXPECTED_CIVIC_CLINICALTRIAL_SIZES);
                LOGGER.warn(keysClinicalTrials);
            }

            clinicalTrialList.add(ImmutableCivicClinicalTrial.builder()
                    .nctId(clinicalTrial.getAsJsonObject().getAsJsonPrimitive("nct_id").getAsString())
                    .description(clinicalTrial.getAsJsonObject().getAsJsonPrimitive("description").getAsString())
                    .clinicalTrialUrl(clinicalTrial.getAsJsonObject().getAsJsonPrimitive("clinical_trial_url").getAsString())
                    .name(clinicalTrial.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .build());
        }
        return clinicalTrialList;
    }

    @NotNull
    private static CivicPublicationDate createPublicationDate(@NotNull JsonObject objectPublicationDate) {
        Set<String> keysPublicationDate = objectPublicationDate.keySet();
        if (!EXPECTED_CIVIC_PUBLICATIONS_DATE_SIZES.contains(keysPublicationDate.size())) {
            LOGGER.warn("Found {} in civic publication date rather than the expected {}",
                    keysPublicationDate.size(),
                    EXPECTED_CIVIC_PUBLICATIONS_DATE_SIZES);
            LOGGER.warn(keysPublicationDate);
        }

        return ImmutableCivicPublicationDate.builder()
                .year(!objectPublicationDate.has("year") ? null : objectPublicationDate.getAsJsonPrimitive("year").getAsString())
                .day(!objectPublicationDate.has("day") ? null : objectPublicationDate.getAsJsonPrimitive("day").getAsString())
                .month(!objectPublicationDate.has("month") ? null : objectPublicationDate.getAsJsonPrimitive("month").getAsString())
                .build();
    }

    @NotNull
    private static CivicDisease createDiseases(@NotNull JsonObject objectDisease) {
        Set<String> keysDisease = objectDisease.keySet();
        if (!EXPECTED_CIVIC_DISEASES_SIZES.contains(keysDisease.size())) {
            LOGGER.warn("Found {} in civic diseases rather than the expected {}", keysDisease.size(), EXPECTED_CIVIC_DISEASES_SIZES);
            LOGGER.warn(keysDisease);
        }

        return ImmutableCivicDisease.builder()
                .doid(objectDisease.get("doid").isJsonNull() ? null : objectDisease.getAsJsonPrimitive("doid").getAsString())
                .url(objectDisease.getAsJsonPrimitive("url").getAsString())
                .displayName(objectDisease.getAsJsonPrimitive("display_name").getAsString())
                .id(objectDisease.getAsJsonPrimitive("id").getAsString())
                .name(objectDisease.getAsJsonPrimitive("name").getAsString())
                .build();
    }

    @NotNull
    private static List<CivicDrug> createDrugs(@NotNull JsonArray arrayDrugs) {
        List<CivicDrug> drugsList = Lists.newArrayList();
        for (JsonElement drug : arrayDrugs) {
            Set<String> keysDrugs = drug.getAsJsonObject().keySet();
            if (!EXPECTED_CIVIC_DRUGS_SIZES.contains(keysDrugs.size())) {
                LOGGER.warn("Found {} in civic drugs rather than the expected {}", keysDrugs.size(), EXPECTED_CIVIC_DRUGS_SIZES);
                LOGGER.warn(keysDrugs);
            }

            drugsList.add(ImmutableCivicDrug.builder()
                    .pubchemId(drug.getAsJsonObject().get("pubchem_id").isJsonNull()
                            ? null
                            : drug.getAsJsonObject().getAsJsonPrimitive("pubchem_id").getAsString())
                    .id(drug.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .name(drug.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .build());
        }
        return drugsList;
    }

    @NotNull
    private static CivicLifecycleActions createLifeCycleActions(@NotNull JsonObject objectLifeCycleActions) {
        Set<String> keysLifecycleActions = objectLifeCycleActions.keySet();
        if (!EXPECTED_CIVIC_LIFECYCLE_ACTIONS_SIZES.contains(keysLifecycleActions.size())) {
            LOGGER.warn("Found {} in civic lifecycle actions rather than the expected {}",
                    keysLifecycleActions.size(),
                    EXPECTED_CIVIC_LIFECYCLE_ACTIONS_SIZES);
            LOGGER.warn(keysLifecycleActions);
        }

        return ImmutableCivicLifecycleActions.builder()
                .lastCommentedOn(objectLifeCycleActions.getAsJsonObject("last_commented_on") == null
                        ? null
                        : createLastCommentedOn(objectLifeCycleActions.getAsJsonObject("last_commented_on")))
                .lastModified(objectLifeCycleActions.getAsJsonObject("last_modified") == null
                        ? null
                        : createLastModified(objectLifeCycleActions.getAsJsonObject("last_modified")))
                .lastReviewed(objectLifeCycleActions.getAsJsonObject("last_reviewed") == null
                        ? null
                        : createLastReviewed(objectLifeCycleActions.getAsJsonObject("last_reviewed")))
                .build();
    }

    @NotNull
    private static CivicLastCommentedOn createLastCommentedOn(@NotNull JsonObject objectLastCommentedOn) {
        Set<String> keysLastCommentedOn = objectLastCommentedOn.keySet();
        if (!EXPECTED_CIVIC_LAST_COMMENTED_ON_SIZES.contains(keysLastCommentedOn.size())) {
            LOGGER.warn("Found {} in civic last commented on rather than the expected {}",
                    keysLastCommentedOn.size(),
                    EXPECTED_CIVIC_LAST_COMMENTED_ON_SIZES);
            LOGGER.warn(keysLastCommentedOn);
        }

        return ImmutableCivicLastCommentedOn.builder()
                .timestamp(objectLastCommentedOn.getAsJsonPrimitive("timestamp").getAsString())
                .user(createCivicUser(objectLastCommentedOn.getAsJsonObject("user")))
                .build();
    }

    @NotNull
    private static CivicLastModified createLastModified(@NotNull JsonObject objectLastModified) {
        Set<String> keysLastModified = objectLastModified.keySet();
        if (!EXPECTED_CIVIC_LAST_MODIFIED_SIZES.contains(keysLastModified.size())) {
            LOGGER.warn("Found {} in civic last modified rather than the expected {}",
                    keysLastModified.size(),
                    EXPECTED_CIVIC_LAST_MODIFIED_SIZES);
            LOGGER.warn(keysLastModified);
        }

        return ImmutableCivicLastModified.builder()
                .timestamp(objectLastModified.getAsJsonPrimitive("timestamp").getAsString())
                .user(createCivicUser(objectLastModified.getAsJsonObject("user")))
                .build();
    }

    @NotNull
    private static CivicLastReviewed createLastReviewed(@NotNull JsonObject objectLastReviewed) {
        Set<String> keysLastReviewed = objectLastReviewed.keySet();
        if (!EXPECTED_CIVIC_LAST_REVIEWED_SIZES.contains(keysLastReviewed.size())) {
            LOGGER.warn("Found {} in civic last reviewed rather than the expected {}",
                    keysLastReviewed.size(),
                    EXPECTED_CIVIC_LAST_REVIEWED_SIZES);
            LOGGER.warn(keysLastReviewed);
        }

        return ImmutableCivicLastReviewed.builder()
                .timestamp(objectLastReviewed.getAsJsonPrimitive("timestamp").getAsString())
                .user(createCivicUser(objectLastReviewed.getAsJsonObject("user")))
                .build();
    }

    @NotNull
    private static CivicUser createCivicUser(@NotNull JsonObject objectUser) {
        Set<String> keysUser = objectUser.keySet();
        if (!EXPECTED_CIVIC_USER_SIZES.contains(keysUser.size())) {
            LOGGER.warn("Found {} in civic user rather than the expected {}", keysUser.size(), EXPECTED_CIVIC_USER_SIZES);
            LOGGER.warn(keysUser);
        }

        return ImmutableCivicUser.builder()
                .username(objectUser.getAsJsonPrimitive("username").getAsString())
                .areaOfExpertise(objectUser.get("area_of_expertise").isJsonNull()
                        ? null
                        : objectUser.getAsJsonPrimitive("area_of_expertise").getAsString())
                .organization(createOrganization(objectUser.getAsJsonObject("organization")))
                .twitterHandle(objectUser.get("twitter_handle").isJsonNull()
                        ? null
                        : objectUser.getAsJsonPrimitive("twitter_handle").getAsString())
                .name(objectUser.getAsJsonPrimitive("name").getAsString())
                .bio(objectUser.get("bio").isJsonNull() ? null : objectUser.getAsJsonPrimitive("bio").getAsString())
                .url(objectUser.get("url").isJsonNull() ? null : objectUser.getAsJsonPrimitive("url").getAsString())
                .createdAt(objectUser.getAsJsonPrimitive("created_at").getAsString())
                .avatars(createAvatars(objectUser.getAsJsonObject("avatars")))
                .acceptedLicense(objectUser.get("accepted_license").isJsonNull()
                        ? null
                        : objectUser.getAsJsonPrimitive("accepted_license").getAsString())
                .affiliation(objectUser.get("affiliation").isJsonNull() ? null : objectUser.getAsJsonPrimitive("affiliation").getAsString())
                .avatarUrl(objectUser.getAsJsonPrimitive("avatar_url").getAsString())
                .role(objectUser.getAsJsonPrimitive("role").getAsString())
                .facebookProfile(objectUser.get("facebook_profile").isJsonNull()
                        ? null
                        : objectUser.getAsJsonPrimitive("facebook_profile").getAsString())
                .linkedinProfile(objectUser.get("linkedin_profile").isJsonNull()
                        ? null
                        : objectUser.getAsJsonPrimitive("linkedin_profile").getAsString())
                .orcid(objectUser.get("orcid").isJsonNull() ? null : objectUser.getAsJsonPrimitive("orcid").getAsString())
                .displayName(objectUser.getAsJsonPrimitive("display_name").getAsString())
                .lastSeenAt(objectUser.get("last_seen_at").isJsonNull()
                        ? null
                        : objectUser.getAsJsonPrimitive("last_seen_at").getAsString())
                .featuredExpert(objectUser.getAsJsonPrimitive("featured_expert").getAsString())
                .id(objectUser.getAsJsonPrimitive("id").getAsString())
                .signupComplete(objectUser.get("signup_complete").isJsonNull()
                        ? null
                        : objectUser.getAsJsonPrimitive("signup_complete").getAsString())
                .build();
    }

    @NotNull
    private static CivicOrganization createOrganization(@NotNull JsonObject objectOrganization) {
        Set<String> keysOrganization = objectOrganization.keySet();
        if (!EXPECTED_CIVIC_ORGANIZATION_SIZES.contains(keysOrganization.size())) {
            LOGGER.warn("Found {} in civic organization rather than the expected {}",
                    keysOrganization.size(),
                    EXPECTED_CIVIC_ORGANIZATION_SIZES);
            LOGGER.warn(keysOrganization);
        }

        return ImmutableCivicOrganization.builder()
                .url(!objectOrganization.has("url") ? null : objectOrganization.getAsJsonPrimitive("url").getAsString())
                .id(!objectOrganization.has("id") ? null : objectOrganization.getAsJsonPrimitive("id").getAsString())
                .profileImage(!objectOrganization.has("profile_image")
                        ? null
                        : createProfileImage(objectOrganization.getAsJsonObject("profile_image")))
                .description(!objectOrganization.has("description")
                        ? null
                        : objectOrganization.getAsJsonPrimitive("description").getAsString())
                .name(!objectOrganization.has("name") ? null : objectOrganization.getAsJsonPrimitive("name").getAsString())
                .build();
    }

    @NotNull
    private static CivicProfileImage createProfileImage(@NotNull JsonObject objectProfileImage) {
        Set<String> keysProfileImage = objectProfileImage.keySet();
        if (!EXPECTED_CIVIC_PROFILE_IMAGE_SIZES.contains(keysProfileImage.size())) {
            LOGGER.warn("Found {} in civic profile image rather than the expected {}",
                    keysProfileImage.size(),
                    EXPECTED_CIVIC_PROFILE_IMAGE_SIZES);
            LOGGER.warn(keysProfileImage);
        }

        return ImmutableCivicProfileImage.builder()
                .x32(objectProfileImage.getAsJsonPrimitive("x32").getAsString())
                .x256(objectProfileImage.getAsJsonPrimitive("x256").getAsString())
                .x14(objectProfileImage.getAsJsonPrimitive("x14").getAsString())
                .x64(objectProfileImage.getAsJsonPrimitive("x64").getAsString())
                .x128(objectProfileImage.getAsJsonPrimitive("x128").getAsString())
                .build();
    }

    @NotNull
    private static CivicAvatars createAvatars(@NotNull JsonObject objectAvatars) {
        Set<String> keysAvatars = objectAvatars.keySet();
        if (!EXPECTED_CIVIC_AVATARS_SIZES.contains(keysAvatars.size())) {
            LOGGER.warn("Found {} in civic avatars rather than the expected {}", keysAvatars.size(), EXPECTED_CIVIC_AVATARS_SIZES);
            LOGGER.warn(keysAvatars);
        }

        return ImmutableCivicAvatars.builder()
                .x32(objectAvatars.getAsJsonPrimitive("x32").getAsString())
                .x14(objectAvatars.getAsJsonPrimitive("x14").getAsString())
                .x64(objectAvatars.getAsJsonPrimitive("x64").getAsString())
                .x128(objectAvatars.getAsJsonPrimitive("x128").getAsString())
                .build();
    }

    @NotNull
    private static List<CivicVariantType> createVariantTypes(@NotNull JsonArray arrayVariantTypes) {
        List<CivicVariantType> civicVariantTypeList = Lists.newArrayList();
        for (JsonElement variantTypes : arrayVariantTypes) {
            Set<String> keysVariantTypes = variantTypes.getAsJsonObject().keySet();
            if (!EXPECTED_CIVIC_VARIANTTYPES_SIZES.contains(keysVariantTypes.size())) {
                LOGGER.warn("Found {} in civic variant types rather than the expected {}",
                        keysVariantTypes.size(),
                        EXPECTED_CIVIC_VARIANTTYPES_SIZES);
                LOGGER.warn(keysVariantTypes);
            }

            civicVariantTypeList.add(ImmutableCivicVariantType.builder()
                    .displayName(variantTypes.getAsJsonObject().getAsJsonPrimitive("display_name").getAsString())
                    .description(variantTypes.getAsJsonObject().getAsJsonPrimitive("description").getAsString())
                    .url(variantTypes.getAsJsonObject().getAsJsonPrimitive("url").getAsString())
                    .soId(variantTypes.getAsJsonObject().getAsJsonPrimitive("so_id").getAsString())
                    .id(variantTypes.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .name(variantTypes.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .build());
        }
        return civicVariantTypeList;
    }
}
