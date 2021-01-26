package com.hartwig.hmftools.ckb.therapy;

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
import com.hartwig.hmftools.ckb.common.ImmutableIndicationInfo;
import com.hartwig.hmftools.ckb.common.ImmutableMolecularProfileInfo;
import com.hartwig.hmftools.ckb.common.ImmutableReferenceInfo;
import com.hartwig.hmftools.ckb.common.ImmutableTherapyInfo;
import com.hartwig.hmftools.ckb.common.IndicationInfo;
import com.hartwig.hmftools.ckb.common.MolecularProfileInfo;
import com.hartwig.hmftools.ckb.common.ReferenceInfo;
import com.hartwig.hmftools.ckb.common.TherapyInfo;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;
import com.hartwig.hmftools.common.utils.json.JsonFunctions;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class TherapyFactory {

    private static final Logger LOGGER = LogManager.getLogger(TherapyFactory.class);

    private TherapyFactory() {

    }

    @NotNull
    public static List<Therapy> readingTherapy(@NotNull String therapyDir) throws IOException {
        LOGGER.info("Start reading therapies");

        List<Therapy> therapies = Lists.newArrayList();
        File[] filesTherapies = new File(therapyDir).listFiles();

        if (filesTherapies != null) {
            LOGGER.info("The total files in the therapy dir is {}", filesTherapies.length);

            for (File therapy : filesTherapies) {
                JsonParser parser = new JsonParser();
                JsonReader reader = new JsonReader(new FileReader(therapy));
                reader.setLenient(true);

                while (reader.peek() != JsonToken.END_DOCUMENT) {
                    JsonObject therapyEntryObject = parser.parse(reader).getAsJsonObject();
                    JsonDatamodelChecker therapyObjectChecker = TherapyDataModelChecker.therapyObjectChecker();
                    therapyObjectChecker.check(therapyEntryObject);

                    therapies.add(ImmutableTherapy.builder()
                            .id(JsonFunctions.string(therapyEntryObject, "id"))
                            .therapyName(JsonFunctions.string(therapyEntryObject, "therapyName"))
                            .synonyms(JsonFunctions.nullableString(therapyEntryObject, "synonyms"))
                            .therapyDescription(createDescription(therapyEntryObject.getAsJsonArray("therapyDescriptions")))
                            .createDate(JsonFunctions.string(therapyEntryObject, "createDate"))
                            .updateDate(JsonFunctions.nullableString(therapyEntryObject, "updateDate"))
                            .evidence(extractEvidence(therapyEntryObject.getAsJsonArray("evidence")))
                            .clinicalTrial(extractClinicalTrial(therapyEntryObject.getAsJsonArray("clinicalTrials")))
                            .drug(extractDrug(therapyEntryObject.getAsJsonArray("drugs")))
                            .globalApprovalStatus(extractGlobalApprovalStatus(therapyEntryObject.getAsJsonArray("globalApprovalStatus")))
                            .build());
                }
                reader.close();
            }
        }
        LOGGER.info("Finished reading therapy");

        return therapies;
    }

    @NotNull
    public static List<TherapyTherapyDescription> createDescription(@NotNull JsonArray jsonArray) {
        List<TherapyTherapyDescription> descriptions = Lists.newArrayList();
        JsonDatamodelChecker descriptionChecker = TherapyDataModelChecker.descriptionObjectChecker();

        for (JsonElement description : jsonArray) {
            JsonObject descriptionJsonObject = description.getAsJsonObject();
            descriptionChecker.check(descriptionJsonObject);

            descriptions.add(ImmutableTherapyTherapyDescription.builder()
                    .description(JsonFunctions.nullableString(descriptionJsonObject, "description"))
                    .reference(extractReferences(descriptionJsonObject.getAsJsonArray("references")))
                    .build());
        }

        return descriptions;
    }

    @NotNull
    public static List<ReferenceInfo> extractReferences(@NotNull JsonArray jsonArray) {
        List<ReferenceInfo> references = Lists.newArrayList();
        JsonDatamodelChecker referenceChecker = TherapyDataModelChecker.referencesObjectChecker();

        for (JsonElement reference : jsonArray) {
            JsonObject referenceJsonObject = reference.getAsJsonObject();
            referenceChecker.check(referenceJsonObject);

            references.add(ImmutableReferenceInfo.builder()
                    .id(JsonFunctions.string(referenceJsonObject, "id"))
                    .pubMedId(JsonFunctions.nullableString(referenceJsonObject, "pubMedId"))
                    .title(JsonFunctions.nullableString(referenceJsonObject, "title"))
                    .url(JsonFunctions.nullableString(referenceJsonObject, "url"))
                    .build());
        }

        return references;
    }

    @NotNull
    public static List<TherapyEvidence> extractEvidence(@NotNull JsonArray jsonArray) {
        List<TherapyEvidence> evidences = Lists.newArrayList();
        JsonDatamodelChecker evidenceChecker = TherapyDataModelChecker.evidenceObjectChecker();

        for (JsonElement evidence : jsonArray) {
            JsonObject evidenceJsonObject = evidence.getAsJsonObject();
            evidenceChecker.check(evidenceJsonObject);

            evidences.add(ImmutableTherapyEvidence.builder()
                    .id(JsonFunctions.string(evidenceJsonObject, "id"))
                    .approvalStatus(JsonFunctions.string(evidenceJsonObject, "approvalStatus"))
                    .evidenceType(JsonFunctions.string(evidenceJsonObject, "evidenceType"))
                    .efficacyEvidence(JsonFunctions.string(evidenceJsonObject, "efficacyEvidence"))
                    .molecularProfile(extractMolecularProfile(evidenceJsonObject.getAsJsonObject("molecularProfile")))
                    .therapy(extractTherapy(evidenceJsonObject.getAsJsonObject("therapy")))
                    .indication(extractIndication(evidenceJsonObject.getAsJsonObject("indication")))
                    .responseType(JsonFunctions.string(evidenceJsonObject, "responseType"))
                    .reference(extractReference(evidenceJsonObject.getAsJsonArray("references")))
                    .ampCapAscoEvidenceLevel(JsonFunctions.string(evidenceJsonObject, "ampCapAscoEvidenceLevel"))
                    .ampCapAscoInferredTier(JsonFunctions.string(evidenceJsonObject, "ampCapAscoInferredTier"))
                    .build());
        }
        return evidences;
    }

    @NotNull
    public static MolecularProfileInfo extractMolecularProfile(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker molecularProfileChecker = TherapyDataModelChecker.molecularProfileObjectChecker();
        molecularProfileChecker.check(jsonObject);

        return ImmutableMolecularProfileInfo.builder()
                .id(JsonFunctions.string(jsonObject, "id"))
                .profileName(JsonFunctions.string(jsonObject, "profileName"))
                .build();
    }

    @NotNull
    public static TherapyInfo extractTherapy(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker therapyChecker = TherapyDataModelChecker.therapyChecker();
        therapyChecker.check(jsonObject);

        return ImmutableTherapyInfo.builder()
                .id(JsonFunctions.string(jsonObject, "id"))
                .therapyName(JsonFunctions.string(jsonObject, "therapyName"))
                .synonyms(JsonFunctions.nullableString(jsonObject, "synonyms"))
                .build();
    }

    @NotNull
    public static IndicationInfo extractIndication(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker indicationChecker = TherapyDataModelChecker.indicationChecker();
        indicationChecker.check(jsonObject);

        return ImmutableIndicationInfo.builder()
                .id(JsonFunctions.string(jsonObject, "id"))
                .name(JsonFunctions.string(jsonObject, "name"))
                .source(JsonFunctions.string(jsonObject, "source"))
                .build();
    }

    @NotNull
    public static List<ReferenceInfo> extractReference(@NotNull JsonArray jsonArray) {
        List<ReferenceInfo> references = Lists.newArrayList();
        JsonDatamodelChecker referenceChecker = TherapyDataModelChecker.referenceChecker();

        for (JsonElement reference : jsonArray) {
            JsonObject referenceJsonObject = reference.getAsJsonObject();
            referenceChecker.check(referenceJsonObject);

            references.add(ImmutableReferenceInfo.builder()
                    .id(JsonFunctions.string(referenceJsonObject, "id"))
                    .pubMedId(JsonFunctions.nullableString(referenceJsonObject, "pubMedId"))
                    .title(JsonFunctions.nullableString(referenceJsonObject, "title"))
                    .url(JsonFunctions.nullableString(referenceJsonObject, "url"))
                    .build());
        }
        return references;
    }

    @NotNull
    public static List<TherapyClinicalTrial> extractClinicalTrial(@NotNull JsonArray jsonArray) {
        List<TherapyClinicalTrial> clinicaltrials = Lists.newArrayList();
        JsonDatamodelChecker clinicalTrialChecker = TherapyDataModelChecker.clinicalTrialChecker();

        for (JsonElement clinicalTrial : jsonArray) {
            JsonObject clinicalTrialJsonObject = clinicalTrial.getAsJsonObject();
            clinicalTrialChecker.check(clinicalTrialJsonObject);

            clinicaltrials.add(ImmutableTherapyClinicalTrial.builder()
                    .nctId(JsonFunctions.string(clinicalTrialJsonObject, "nctId"))
                    .title(JsonFunctions.string(clinicalTrialJsonObject, "title"))
                    .phase(JsonFunctions.string(clinicalTrialJsonObject, "phase"))
                    .recruitment(JsonFunctions.string(clinicalTrialJsonObject, "recruitment"))
                    .therapy(extractTherapyList(clinicalTrialJsonObject.getAsJsonArray("therapies")))
                    .build());
        }
        return clinicaltrials;
    }

    @NotNull
    public static List<TherapyInfo> extractTherapyList(@NotNull JsonArray jsonArray) {
        List<TherapyInfo> therapies = Lists.newArrayList();
        JsonDatamodelChecker therapyChecker = TherapyDataModelChecker.therapyChecker();

        for (JsonElement therapy : jsonArray) {
            JsonObject therapyJsonObject = therapy.getAsJsonObject();
            therapyChecker.check(therapyJsonObject);

            therapies.add(ImmutableTherapyInfo.builder()
                    .id(JsonFunctions.string(therapyJsonObject, "id"))
                    .therapyName(JsonFunctions.string(therapyJsonObject, "therapyName"))
                    .synonyms(JsonFunctions.nullableString(therapyJsonObject, "synonyms"))
                    .build());
        }
        return therapies;
    }

    @NotNull
    public static List<TherapyDrug> extractDrug(@NotNull JsonArray jsonArray) {
        List<TherapyDrug> drugs = Lists.newArrayList();
        JsonDatamodelChecker drugChecker = TherapyDataModelChecker.drugsChecker();

        for (JsonElement drug : jsonArray) {
            JsonObject drugJsonObject = drug.getAsJsonObject();
            drugChecker.check(drugJsonObject);

            drugs.add(ImmutableTherapyDrug.builder()
                    .id(JsonFunctions.string(drugJsonObject, "id"))
                    .drugName(JsonFunctions.string(drugJsonObject, "drugName"))
                    .terms(JsonFunctions.stringList(drugJsonObject, "terms"))
                    .build());

        }
        return drugs;
    }

    @NotNull
    public static List<TherapyGlobalApprovalStatus> extractGlobalApprovalStatus(@NotNull JsonArray jsonArray) {
        List<TherapyGlobalApprovalStatus> globalApprovalStatuses = Lists.newArrayList();
        JsonDatamodelChecker globalApprovalStatusChecker = TherapyDataModelChecker.globalApprovalStatusChecker();

        for (JsonElement globalApprovalStatus : jsonArray) {
            JsonObject globalApprovalStatusJsonObject = globalApprovalStatus.getAsJsonObject();
            globalApprovalStatusChecker.check(globalApprovalStatusJsonObject);

            globalApprovalStatuses.add(ImmutableTherapyGlobalApprovalStatus.builder()
                    .id(JsonFunctions.string(globalApprovalStatusJsonObject, "id"))
                    .therapy(extractTherapy(globalApprovalStatusJsonObject.getAsJsonObject("therapy")))
                    .indication(extractIndication(globalApprovalStatusJsonObject.getAsJsonObject("indication")))
                    .molecularProfile(extractMolecularProfile(globalApprovalStatusJsonObject.getAsJsonObject("molecularProfile")))
                    .approvalAuthority(JsonFunctions.string(globalApprovalStatusJsonObject, "approvalAuthority"))
                    .approvalStatus(JsonFunctions.string(globalApprovalStatusJsonObject, "approvalStatus"))
                    .build());
        }
        return globalApprovalStatuses;

    }
}
