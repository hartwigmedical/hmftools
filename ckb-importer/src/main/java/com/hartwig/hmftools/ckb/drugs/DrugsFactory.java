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
            drugDescriptions.add(ImmutableDrugDescription.builder()
                    .description(JsonFunctions.string(drugDescription.getAsJsonObject(), "description"))
                    .references(extractDrugReferences(drugDescription.getAsJsonObject().getAsJsonArray("references")))
                    .build());
        }
        return drugDescriptions;
    }

    @NotNull
    public static List<DrugsReferences> extractDrugReferences(@NotNull JsonArray jsonArray) {
        List<DrugsReferences> drugReferences = Lists.newArrayList();
        for (JsonElement drugReference : jsonArray) {
            drugReferences.add(ImmutableDrugsReferences.builder()
                    .id(JsonFunctions.string(drugReference.getAsJsonObject(), "id"))
                    .pubMedId(JsonFunctions.nullableString(drugReference.getAsJsonObject(), "pubMedId"))
                    .title(JsonFunctions.nullableString(drugReference.getAsJsonObject(), "title"))
                    .url(JsonFunctions.nullableString(drugReference.getAsJsonObject(), "url"))
                    .authors(JsonFunctions.nullableString(drugReference.getAsJsonObject(), "authors"))
                    .journal(JsonFunctions.nullableString(drugReference.getAsJsonObject(), "journal"))
                    .volume(JsonFunctions.nullableString(drugReference.getAsJsonObject(), "volume"))
                    .issue(JsonFunctions.nullableString(drugReference.getAsJsonObject(), "issue"))
                    .date(JsonFunctions.nullableString(drugReference.getAsJsonObject(), "date"))
                    .abstractText(JsonFunctions.nullableString(drugReference.getAsJsonObject(), "abstractText"))
                    .year(JsonFunctions.nullableString(drugReference.getAsJsonObject(), "year"))
                    .build());
        }
        return drugReferences;
    }

    @NotNull
    public static List<DrugClasses> extractDrugsClasses(@NotNull JsonArray jsonArray) {
        List<DrugClasses> drugClasses = Lists.newArrayList();
        for (JsonElement drugClass : jsonArray) {
            drugClasses.add(ImmutableDrugClasses.builder()
                    .id(JsonFunctions.string(drugClass.getAsJsonObject(), "id"))
                    .drugClass(JsonFunctions.string(drugClass.getAsJsonObject(), "drugClass"))
                    .build());
        }
        return drugClasses;
    }

    @NotNull
    public static List<DrugsClinicalTrials> extractCliniclaTrials(@NotNull JsonArray jsonArray) {
        List<DrugsClinicalTrials> clinicalTrials = Lists.newArrayList();
        for (JsonElement clinicalTrial : jsonArray) {
            clinicalTrials.add(ImmutableDrugsClinicalTrials.builder()
                    .nctId(JsonFunctions.string(clinicalTrial.getAsJsonObject(), "nctId"))
                    .title(JsonFunctions.string(clinicalTrial.getAsJsonObject(), "title"))
                    .phase(JsonFunctions.string(clinicalTrial.getAsJsonObject(), "phase"))
                    .recruitment(JsonFunctions.string(clinicalTrial.getAsJsonObject(), "recruitment"))
                    .therapies(extractDrugTherapies(clinicalTrial.getAsJsonObject().getAsJsonArray("therapies")))
                    .build());
        }
        return clinicalTrials;
    }

    @NotNull
    public static List<DrugsTherapies> extractDrugTherapies(@NotNull JsonArray jsonArray) {
        List<DrugsTherapies> drugTherapies = Lists.newArrayList();
        for (JsonElement drugTherpy : jsonArray) {
            drugTherapies.add(ImmutableDrugsTherapies.builder()
                    .id(JsonFunctions.string(drugTherpy.getAsJsonObject(), "id"))
                    .therapyName(JsonFunctions.string(drugTherpy.getAsJsonObject(), "therapyName"))
                    .synonyms(JsonFunctions.nullableString(drugTherpy.getAsJsonObject(), "synonyms"))
                    .build());
        }
        return drugTherapies;
    }

    @NotNull
    public static List<DrugsEvidence> extractEvidence(@NotNull JsonArray jsonArray) {
        List<DrugsEvidence> evidences = Lists.newArrayList();
        for (JsonElement evidence : jsonArray) {
            evidences.add(ImmutableDrugsEvidence.builder()
                    .id(JsonFunctions.string(evidence.getAsJsonObject(), "id"))
                    .approvalStatus(JsonFunctions.string(evidence.getAsJsonObject(), "approvalStatus"))
                    .evidenceType(JsonFunctions.string(evidence.getAsJsonObject(), "evidenceType"))
                    .efficacyEvidence(JsonFunctions.string(evidence.getAsJsonObject(), "efficacyEvidence"))
                    .molecularProfile(extractMolecularProfile(evidence.getAsJsonObject().getAsJsonObject("molecularProfile")))
                    .therapy(extractTherapy(evidence.getAsJsonObject().getAsJsonObject("therapy")))
                    .indications(evidence.getAsJsonObject().has("indications") ? extractIndications(evidence.getAsJsonObject().getAsJsonObject("indications")) : null)
                    .responseType(JsonFunctions.string(evidence.getAsJsonObject(), "responseType"))
                    .references(extractEvidenceReferences(evidence.getAsJsonObject().getAsJsonArray("references")))
                    .ampCapAscoEvidenceLevel(JsonFunctions.string(evidence.getAsJsonObject(), "ampCapAscoEvidenceLevel"))
                    .ampCapAscoInferredTier(JsonFunctions.string(evidence.getAsJsonObject(), "ampCapAscoInferredTier"))
                    .build());
        }
        return evidences;
    }

    @NotNull
    public static DrugsMolecularProfile extractMolecularProfile(@NotNull JsonObject jsonObject) {
        return ImmutableDrugsMolecularProfile.builder()
                .id(JsonFunctions.string(jsonObject.getAsJsonObject(), "id"))
                .profileName(JsonFunctions.string(jsonObject.getAsJsonObject(), "profileName"))
                .build();
    }

    @NotNull
    public static DrugsTherapy extractTherapy(@NotNull JsonObject jsonObject) {
        return ImmutableDrugsTherapy.builder()
                .id(JsonFunctions.string(jsonObject.getAsJsonObject(), "id"))
                .therapyName(JsonFunctions.string(jsonObject.getAsJsonObject(), "therapyName"))
                .synonyms(JsonFunctions.nullableString(jsonObject.getAsJsonObject(), "synonyms"))
                .build();
    }

    @NotNull
    public static DrugsIndications extractIndications(@NotNull JsonObject jsonObject) {
        return ImmutableDrugsIndications.builder()
                .id(JsonFunctions.string(jsonObject.getAsJsonObject(), "id"))
                .name(JsonFunctions.string(jsonObject.getAsJsonObject(), "name"))
                .source(JsonFunctions.string(jsonObject.getAsJsonObject(), "source"))
                .build();
    }

    @NotNull
    public static List<DrugsEvidenceReferences> extractEvidenceReferences(@NotNull JsonArray jsonArray) {
        List<DrugsEvidenceReferences> references = Lists.newArrayList();
        for (JsonElement reference : jsonArray) {
            references.add(ImmutableDrugsEvidenceReferences.builder()
                    .id(JsonFunctions.string(reference.getAsJsonObject(), "id"))
                    .pubMedId(JsonFunctions.nullableString(reference.getAsJsonObject(), "pubMedId"))
                    .title(JsonFunctions.nullableString(reference.getAsJsonObject(), "title"))
                    .url(JsonFunctions.nullableString(reference.getAsJsonObject(), "url"))
                    .build());
        }
        return references;
    }

    @NotNull
    public static List<DrugsTherapies> extractTherapies(@NotNull JsonArray jsonArray) {
        List<DrugsTherapies> drugsTherapies = Lists.newArrayList();
        for (JsonElement drugTherapy : jsonArray) {
            drugsTherapies.add(ImmutableDrugsTherapies.builder()
                    .id(JsonFunctions.string(drugTherapy.getAsJsonObject(), "id"))
                    .therapyName(JsonFunctions.string(drugTherapy.getAsJsonObject(), "therapyName"))
                    .synonyms(JsonFunctions.nullableString(drugTherapy.getAsJsonObject(), "synonyms"))
                    .build());
        }
        return drugsTherapies;
    }

    @NotNull
    public static List<DrugsGlobalApproavalStatus> extractGlobalApprovaStatus(@NotNull JsonArray jsonArray) {
        List<DrugsGlobalApproavalStatus> globalApproavalStatuses = Lists.newArrayList();
        for (JsonElement globalApproavalStatus : jsonArray) {
            globalApproavalStatuses.add(ImmutableDrugsGlobalApproavalStatus.builder()
                    .id(JsonFunctions.string(globalApproavalStatus.getAsJsonObject(), "id"))
                    .therapy(extractTherapiesGlobalApprovalStatus(globalApproavalStatus.getAsJsonObject().getAsJsonObject("therapy")))
                    .indications(extractIndicationsGlobalApprovalStatus(globalApproavalStatus.getAsJsonObject().getAsJsonObject("indications")))
                    .molecularProfile(extractMolecularProfileGlobalApprovalStatus(globalApproavalStatus.getAsJsonObject().getAsJsonObject("molecularProfile")))
                    .approvalAuthority(JsonFunctions.string(globalApproavalStatus.getAsJsonObject(), "approvalAuthority"))
                    .approvalStatus(JsonFunctions.string(globalApproavalStatus.getAsJsonObject(), "approvalStatus"))
                    .build());
        }
        return globalApproavalStatuses;
    }

    @NotNull
    public static DrugsTherapy extractTherapiesGlobalApprovalStatus(@NotNull JsonObject jsonObject) {
        return ImmutableDrugsTherapy.builder()
                .id(JsonFunctions.string(jsonObject.getAsJsonObject(), "id"))
                .therapyName(JsonFunctions.string(jsonObject.getAsJsonObject(), "therapyName"))
                .synonyms(JsonFunctions.string(jsonObject.getAsJsonObject(), "synonyms"))
                .build();
    }

    @NotNull
    public static DrugsIndications extractIndicationsGlobalApprovalStatus(@NotNull JsonObject jsonObject) {
        return ImmutableDrugsIndications.builder()
                .id(JsonFunctions.string(jsonObject.getAsJsonObject(), "id"))
                .name(JsonFunctions.string(jsonObject.getAsJsonObject(), "name"))
                .source(JsonFunctions.string(jsonObject.getAsJsonObject(), "source"))
                .build();
    }

    @NotNull
    public static DrugsMolecularProfile extractMolecularProfileGlobalApprovalStatus(@NotNull JsonObject jsonObject) {
        return ImmutableDrugsMolecularProfile.builder()
                .id(JsonFunctions.string(jsonObject.getAsJsonObject(), "id"))
                .profileName(JsonFunctions.string(jsonObject.getAsJsonObject(), "profileName"))
                .build();
    }
}
