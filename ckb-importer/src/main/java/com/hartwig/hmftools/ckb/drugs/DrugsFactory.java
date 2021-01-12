package com.hartwig.hmftools.ckb.drugs;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;
import com.hartwig.hmftools.common.utils.json.JsonFunctions;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class DrugsFactory {

    private static final Logger LOGGER = LogManager.getLogger(DrugsFactory.class);

    private DrugsFactory() {
    }

    @NotNull
    public static List<Drugs> readingDrugs(@NotNull String drugsDir) throws IOException {
        LOGGER.info("Start reading drugs");

        List<Drugs> drugs = Lists.newArrayList();
        File[] filesDrugs = new File(drugsDir).listFiles();

        if (filesDrugs != null) {
            LOGGER.info("The total files in the drugs dir is {}", filesDrugs.length);

            for (File drug : filesDrugs) {
                JsonParser parser = new JsonParser();
                JsonReader reader = new JsonReader(new FileReader(drug));
                reader.setLenient(true);

                while (reader.peek() != JsonToken.END_DOCUMENT) {
                    JsonObject drugsEntryObject = parser.parse(reader).getAsJsonObject();
                    drugs.add(ImmutableDrugs.builder()
                            .id(JsonFunctions.string(drugsEntryObject, "id"))
                            .drugName(JsonFunctions.string(drugsEntryObject, "drugName"))
                            .terms(JsonFunctions.stringList(drugsEntryObject, "terms"))
                            .synonyms(JsonFunctions.stringList(drugsEntryObject, "synonyms"))
                            .tradeName(JsonFunctions.nullableString(drugsEntryObject, "tradeName"))
                            .drugDescriptions(extractDrugDescriptions(drugsEntryObject.getAsJsonArray("drugDescriptions")))
                            .drugClasses(extractDrugsClasses(drugsEntryObject.getAsJsonArray("drugClasses")))
                            .casRegistryNum(JsonFunctions.nullableString(drugsEntryObject, "casRegistryNum"))
                            .nctId(JsonFunctions.nullableString(drugsEntryObject, "ncitId"))
                            .createDate(JsonFunctions.string(drugsEntryObject, "createDate"))
                            .clinicalTrials(extractCliniclaTrials(drugsEntryObject.getAsJsonArray("clinicalTrials")))
                            .evidence(extractEvidence(drugsEntryObject.getAsJsonArray("evidence")))
                            .therapies(extractTherapies(drugsEntryObject.getAsJsonArray("therapies")))
                            .globalApprovaStatus(drugsEntryObject.has("globalApprovaStatus") ? extractGlobalApprovaStatus(drugsEntryObject.getAsJsonArray("globalApprovaStatus")) : null)
                            .build());
                }
            }
        }
        return drugs;
    }

    @NotNull
    public static List<DrugDescription> extractDrugDescriptions(@NotNull JsonArray jsonArray) {
        List<DrugDescription> drugDescriptions = Lists.newArrayList();

        for (JsonElement drugDescription : jsonArray) {
            JsonObject drugDescriptionObject = drugDescription.getAsJsonObject();
            drugDescriptions.add(ImmutableDrugDescription.builder()
                    .description(JsonFunctions.string(drugDescriptionObject, "description"))
                    .references(extractDrugReferences(drugDescriptionObject.getAsJsonArray("references")))
                    .build());
        }
        return drugDescriptions;
    }

    @NotNull
    public static List<DrugsReferences> extractDrugReferences(@NotNull JsonArray jsonArray) {
        List<DrugsReferences> drugReferences = Lists.newArrayList();
        for (JsonElement drugReference : jsonArray) {
            JsonObject drugsReferenceObject = drugReference.getAsJsonObject();
            drugReferences.add(ImmutableDrugsReferences.builder()
                    .id(JsonFunctions.string(drugsReferenceObject, "id"))
                    .pubMedId(JsonFunctions.nullableString(drugsReferenceObject, "pubMedId"))
                    .title(JsonFunctions.nullableString(drugsReferenceObject, "title"))
                    .url(JsonFunctions.nullableString(drugsReferenceObject, "url"))
                    .authors(JsonFunctions.nullableString(drugsReferenceObject, "authors"))
                    .journal(JsonFunctions.nullableString(drugsReferenceObject, "journal"))
                    .volume(JsonFunctions.nullableString(drugsReferenceObject, "volume"))
                    .issue(JsonFunctions.nullableString(drugsReferenceObject, "issue"))
                    .date(JsonFunctions.nullableString(drugsReferenceObject, "date"))
                    .abstractText(JsonFunctions.nullableString(drugsReferenceObject, "abstractText"))
                    .year(JsonFunctions.nullableString(drugsReferenceObject, "year"))
                    .build());
        }
        return drugReferences;
    }

    @NotNull
    public static List<DrugClasses> extractDrugsClasses(@NotNull JsonArray jsonArray) {
        List<DrugClasses> drugClasses = Lists.newArrayList();
        for (JsonElement drugClass : jsonArray) {
            JsonObject drugClassObject = drugClass.getAsJsonObject();
            drugClasses.add(ImmutableDrugClasses.builder()
                    .id(JsonFunctions.string(drugClassObject, "id"))
                    .drugClass(JsonFunctions.string(drugClassObject, "drugClass"))
                    .build());
        }
        return drugClasses;
    }

    @NotNull
    public static List<DrugsClinicalTrials> extractCliniclaTrials(@NotNull JsonArray jsonArray) {
        List<DrugsClinicalTrials> clinicalTrials = Lists.newArrayList();
        for (JsonElement clinicalTrial : jsonArray) {
            JsonObject clinicalTrialObject = clinicalTrial.getAsJsonObject();
            clinicalTrials.add(ImmutableDrugsClinicalTrials.builder()
                    .nctId(JsonFunctions.string(clinicalTrialObject, "nctId"))
                    .title(JsonFunctions.string(clinicalTrialObject, "title"))
                    .phase(JsonFunctions.string(clinicalTrialObject, "phase"))
                    .recruitment(JsonFunctions.string(clinicalTrialObject, "recruitment"))
                    .therapies(extractDrugTherapies(clinicalTrialObject.getAsJsonArray("therapies")))
                    .build());
        }
        return clinicalTrials;
    }

    @NotNull
    public static List<DrugsTherapies> extractDrugTherapies(@NotNull JsonArray jsonArray) {
        List<DrugsTherapies> drugTherapies = Lists.newArrayList();
        for (JsonElement drugTherpy : jsonArray) {
            JsonObject drugTherapyObject = drugTherpy.getAsJsonObject();
            drugTherapies.add(ImmutableDrugsTherapies.builder()
                    .id(JsonFunctions.string(drugTherapyObject, "id"))
                    .therapyName(JsonFunctions.string(drugTherapyObject, "therapyName"))
                    .synonyms(JsonFunctions.nullableString(drugTherapyObject, "synonyms"))
                    .build());
        }
        return drugTherapies;
    }

    @NotNull
    public static List<DrugsEvidence> extractEvidence(@NotNull JsonArray jsonArray) {
        List<DrugsEvidence> evidences = Lists.newArrayList();
        for (JsonElement evidence : jsonArray) {
            JsonObject evidenceObject = evidence.getAsJsonObject();
            evidences.add(ImmutableDrugsEvidence.builder()
                    .id(JsonFunctions.string(evidenceObject, "id"))
                    .approvalStatus(JsonFunctions.string(evidenceObject, "approvalStatus"))
                    .evidenceType(JsonFunctions.string(evidenceObject, "evidenceType"))
                    .efficacyEvidence(JsonFunctions.string(evidenceObject, "efficacyEvidence"))
                    .molecularProfile(extractMolecularProfile(evidenceObject.getAsJsonObject("molecularProfile")))
                    .therapy(extractTherapy(evidenceObject.getAsJsonObject("therapy")))
                    .indications(evidenceObject.has("indications") ? extractIndications(evidenceObject.getAsJsonObject("indications")) : null)
                    .responseType(JsonFunctions.string(evidenceObject, "responseType"))
                    .references(extractEvidenceReferences(evidenceObject.getAsJsonArray("references")))
                    .ampCapAscoEvidenceLevel(JsonFunctions.string(evidenceObject, "ampCapAscoEvidenceLevel"))
                    .ampCapAscoInferredTier(JsonFunctions.string(evidenceObject, "ampCapAscoInferredTier"))
                    .build());
        }
        return evidences;
    }

    @NotNull
    public static DrugsMolecularProfile extractMolecularProfile(@NotNull JsonObject jsonObject) {
        JsonObject molecularProfileObject = jsonObject.getAsJsonObject();
        return ImmutableDrugsMolecularProfile.builder()
                .id(JsonFunctions.string(molecularProfileObject, "id"))
                .profileName(JsonFunctions.string(molecularProfileObject, "profileName"))
                .build();
    }

    @NotNull
    public static DrugsTherapy extractTherapy(@NotNull JsonObject jsonObject) {
        JsonObject therapyObject = jsonObject.getAsJsonObject();
        return ImmutableDrugsTherapy.builder()
                .id(JsonFunctions.string(therapyObject, "id"))
                .therapyName(JsonFunctions.string(therapyObject, "therapyName"))
                .synonyms(JsonFunctions.nullableString(therapyObject, "synonyms"))
                .build();
    }

    @NotNull
    public static DrugsIndications extractIndications(@NotNull JsonObject jsonObject) {
        JsonObject indicationObject = jsonObject.getAsJsonObject();
        return ImmutableDrugsIndications.builder()
                .id(JsonFunctions.string(indicationObject, "id"))
                .name(JsonFunctions.string(indicationObject, "name"))
                .source(JsonFunctions.string(indicationObject, "source"))
                .build();
    }

    @NotNull
    public static List<DrugsEvidenceReferences> extractEvidenceReferences(@NotNull JsonArray jsonArray) {
        List<DrugsEvidenceReferences> references = Lists.newArrayList();
        for (JsonElement reference : jsonArray) {
            JsonObject referenceObject = reference.getAsJsonObject();
            references.add(ImmutableDrugsEvidenceReferences.builder()
                    .id(JsonFunctions.string(referenceObject, "id"))
                    .pubMedId(JsonFunctions.nullableString(referenceObject, "pubMedId"))
                    .title(JsonFunctions.nullableString(referenceObject, "title"))
                    .url(JsonFunctions.nullableString(referenceObject, "url"))
                    .build());
        }
        return references;
    }

    @NotNull
    public static List<DrugsTherapies> extractTherapies(@NotNull JsonArray jsonArray) {
        List<DrugsTherapies> drugsTherapies = Lists.newArrayList();
        for (JsonElement drugTherapy : jsonArray) {
            JsonObject drugTherapyObject = drugTherapy.getAsJsonObject();
            drugsTherapies.add(ImmutableDrugsTherapies.builder()
                    .id(JsonFunctions.string(drugTherapyObject, "id"))
                    .therapyName(JsonFunctions.string(drugTherapyObject, "therapyName"))
                    .synonyms(JsonFunctions.nullableString(drugTherapyObject, "synonyms"))
                    .build());
        }
        return drugsTherapies;
    }

    @NotNull
    public static List<DrugsGlobalApproavalStatus> extractGlobalApprovaStatus(@NotNull JsonArray jsonArray) {
        List<DrugsGlobalApproavalStatus> globalApproavalStatuses = Lists.newArrayList();
        for (JsonElement globalApproavalStatus : jsonArray) {
            JsonObject globalTherapyObject = globalApproavalStatus.getAsJsonObject();
            globalApproavalStatuses.add(ImmutableDrugsGlobalApproavalStatus.builder()
                    .id(JsonFunctions.string(globalTherapyObject, "id"))
                    .therapy(extractTherapiesGlobalApprovalStatus(globalTherapyObject.getAsJsonObject("therapy")))
                    .indications(extractIndicationsGlobalApprovalStatus(globalTherapyObject.getAsJsonObject("indications")))
                    .molecularProfile(extractMolecularProfileGlobalApprovalStatus(globalTherapyObject.getAsJsonObject("molecularProfile")))
                    .approvalAuthority(JsonFunctions.string(globalTherapyObject, "approvalAuthority"))
                    .approvalStatus(JsonFunctions.string(globalTherapyObject, "approvalStatus"))
                    .build());
        }
        return globalApproavalStatuses;
    }

    @NotNull
    public static DrugsTherapy extractTherapiesGlobalApprovalStatus(@NotNull JsonObject jsonObject) {
        JsonObject therapiesGlobalApprovalStatusObject = jsonObject.getAsJsonObject();
        return ImmutableDrugsTherapy.builder()
                .id(JsonFunctions.string(therapiesGlobalApprovalStatusObject, "id"))
                .therapyName(JsonFunctions.string(therapiesGlobalApprovalStatusObject, "therapyName"))
                .synonyms(JsonFunctions.string(therapiesGlobalApprovalStatusObject, "synonyms"))
                .build();
    }

    @NotNull
    public static DrugsIndications extractIndicationsGlobalApprovalStatus(@NotNull JsonObject jsonObject) {
        JsonObject indicationsGlobalApprovalStatusObject = jsonObject.getAsJsonObject();
        return ImmutableDrugsIndications.builder()
                .id(JsonFunctions.string(indicationsGlobalApprovalStatusObject, "id"))
                .name(JsonFunctions.string(indicationsGlobalApprovalStatusObject, "name"))
                .source(JsonFunctions.string(indicationsGlobalApprovalStatusObject, "source"))
                .build();
    }

    @NotNull
    public static DrugsMolecularProfile extractMolecularProfileGlobalApprovalStatus(@NotNull JsonObject jsonObject) {
        JsonObject molecularProfileGlobalApprovalStatus = jsonObject.getAsJsonObject();
        return ImmutableDrugsMolecularProfile.builder()
                .id(JsonFunctions.string(molecularProfileGlobalApprovalStatus, "id"))
                .profileName(JsonFunctions.string(molecularProfileGlobalApprovalStatus, "profileName"))
                .build();
    }
}
