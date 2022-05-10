package com.hartwig.hmftools.ckb.json.therapy;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.ckb.json.CkbJsonDirectoryReader;
import com.hartwig.hmftools.ckb.json.common.ClinicalTrialInfo;
import com.hartwig.hmftools.ckb.json.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.json.common.DrugInfo;
import com.hartwig.hmftools.ckb.json.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.json.common.GlobalApprovalStatusInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableClinicalTrialInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableDescriptionInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableDrugInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableEvidenceInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableGlobalApprovalStatusInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableIndicationInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableMolecularProfileInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableReferenceInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableTherapyInfo;
import com.hartwig.hmftools.ckb.json.common.IndicationInfo;
import com.hartwig.hmftools.ckb.json.common.MolecularProfileInfo;
import com.hartwig.hmftools.ckb.json.common.ReferenceInfo;
import com.hartwig.hmftools.ckb.json.common.TherapyInfo;
import com.hartwig.hmftools.ckb.util.DateConverter;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;
import com.hartwig.hmftools.common.utils.json.JsonFunctions;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class TherapyReader extends CkbJsonDirectoryReader<JsonTherapy> {
    private static final Logger LOGGER = LogManager.getLogger(TherapyReader.class);

    public TherapyReader(@Nullable final Integer maxFilesToRead) {
        super(maxFilesToRead);
    }

    @NotNull
    @Override
    protected JsonTherapy read(@NotNull final JsonObject object) {
        JsonDatamodelChecker therapyObjectChecker = TherapyDataModelChecker.therapyObjectChecker();
        therapyObjectChecker.check(object);

        return ImmutableJsonTherapy.builder()
                .id(JsonFunctions.integer(object, "id"))
                .therapyName(JsonFunctions.string(object, "therapyName"))
                .synonyms(JsonFunctions.optionalStringList(object, "synonyms"))
                .descriptions(extractDescriptions(object.getAsJsonArray("therapyDescriptions")))
                .createDate(DateConverter.toDate(JsonFunctions.string(object, "createDate")))
                .updateDate(DateConverter.toDate(JsonFunctions.nullableString(object, "updateDate")))
                .evidence(extractEvidence(object.getAsJsonArray("evidence")))
                .clinicalTrials(extractClinicalTrials(object.getAsJsonArray("clinicalTrials")))
                .drugs(extractDrugs(object.getAsJsonArray("drugs")))
                .globalApprovalStatuses(extractGlobalApprovalStatuses(object.getAsJsonArray("globalApprovalStatus")))
                .build();
    }

    @NotNull
    private static List<DescriptionInfo> extractDescriptions(@NotNull JsonArray jsonArray) {
        List<DescriptionInfo> descriptions = Lists.newArrayList();
        JsonDatamodelChecker descriptionChecker = TherapyDataModelChecker.descriptionObjectChecker();

        for (JsonElement description : jsonArray) {
            JsonObject descriptionJsonObject = description.getAsJsonObject();
            descriptionChecker.check(descriptionJsonObject);

            descriptions.add(ImmutableDescriptionInfo.builder()
                    .description(JsonFunctions.nullableString(descriptionJsonObject, "description"))
                    .references(extractReferences(descriptionJsonObject.getAsJsonArray("references")))
                    .build());
        }
        return descriptions;
    }

    @NotNull
    private static List<ReferenceInfo> extractReferences(@NotNull JsonArray jsonArray) {
        List<ReferenceInfo> references = Lists.newArrayList();
        JsonDatamodelChecker referenceChecker = TherapyDataModelChecker.referenceObjectChecker();

        for (JsonElement reference : jsonArray) {
            JsonObject referenceJsonObject = reference.getAsJsonObject();
            referenceChecker.check(referenceJsonObject);

            references.add(ImmutableReferenceInfo.builder()
                    .id(JsonFunctions.integer(referenceJsonObject, "id"))
                    .pubMedId(JsonFunctions.nullableString(referenceJsonObject, "pubMedId"))
                    .title(JsonFunctions.nullableString(referenceJsonObject, "title"))
                    .url(JsonFunctions.nullableString(referenceJsonObject, "url"))
                    .build());
        }
        return references;
    }

    @NotNull
    private static List<EvidenceInfo> extractEvidence(@NotNull JsonArray jsonArray) {
        List<EvidenceInfo> evidences = Lists.newArrayList();
        JsonDatamodelChecker evidenceChecker = TherapyDataModelChecker.evidenceObjectChecker();

        for (JsonElement evidence : jsonArray) {
            JsonObject evidenceJsonObject = evidence.getAsJsonObject();
            evidenceChecker.check(evidenceJsonObject);

            evidences.add(ImmutableEvidenceInfo.builder()
                    .id(JsonFunctions.integer(evidenceJsonObject, "id"))
                    .approvalStatus(JsonFunctions.string(evidenceJsonObject, "approvalStatus"))
                    .evidenceType(JsonFunctions.string(evidenceJsonObject, "evidenceType"))
                    .efficacyEvidence(JsonFunctions.string(evidenceJsonObject, "efficacyEvidence"))
                    .molecularProfile(extractMolecularProfile(evidenceJsonObject.getAsJsonObject("molecularProfile")))
                    .therapy(extractTherapy(evidenceJsonObject.getAsJsonObject("therapy")))
                    .indication(extractIndication(evidenceJsonObject.getAsJsonObject("indication")))
                    .responseType(JsonFunctions.string(evidenceJsonObject, "responseType"))
                    .references(extractReferences(evidenceJsonObject.getAsJsonArray("references")))
                    .ampCapAscoEvidenceLevel(JsonFunctions.string(evidenceJsonObject, "ampCapAscoEvidenceLevel"))
                    .ampCapAscoInferredTier(JsonFunctions.string(evidenceJsonObject, "ampCapAscoInferredTier"))
                    .build());
        }
        return evidences;
    }

    @NotNull
    private static MolecularProfileInfo extractMolecularProfile(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker molecularProfileChecker = TherapyDataModelChecker.molecularProfileObjectChecker();
        molecularProfileChecker.check(jsonObject);

        return ImmutableMolecularProfileInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .profileName(JsonFunctions.string(jsonObject, "profileName"))
                .build();
    }

    @NotNull
    private static TherapyInfo extractTherapy(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker therapyChecker = TherapyDataModelChecker.therapyChecker();
        therapyChecker.check(jsonObject);

        return ImmutableTherapyInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .therapyName(JsonFunctions.string(jsonObject, "therapyName"))
                .synonyms(JsonFunctions.optionalStringList(jsonObject, "synonyms"))
                .build();
    }

    @NotNull
    private static IndicationInfo extractIndication(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker indicationChecker = TherapyDataModelChecker.indicationChecker();
        indicationChecker.check(jsonObject);

        return ImmutableIndicationInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .name(JsonFunctions.string(jsonObject, "name"))
                .source(JsonFunctions.string(jsonObject, "source"))
                .build();
    }

    @NotNull
    private static List<ClinicalTrialInfo> extractClinicalTrials(@NotNull JsonArray jsonArray) {
        List<ClinicalTrialInfo> clinicalTrials = Lists.newArrayList();
        JsonDatamodelChecker clinicalTrialChecker = TherapyDataModelChecker.clinicalTrialChecker();

        for (JsonElement clinicalTrial : jsonArray) {
            JsonObject clinicalTrialJsonObject = clinicalTrial.getAsJsonObject();
            clinicalTrialChecker.check(clinicalTrialJsonObject);

            String nctId = JsonFunctions.string(clinicalTrialJsonObject, "nctId");
            String phase = JsonFunctions.nullableString(clinicalTrialJsonObject, "phase");

            if (phase == null) {
                LOGGER.warn("phase of study '{}' is null in TherapyReader", nctId);
            }

            clinicalTrials.add(ImmutableClinicalTrialInfo.builder()
                    .nctId(nctId)
                    .title(JsonFunctions.string(clinicalTrialJsonObject, "title"))
                    .phase(phase)
                    .recruitment(JsonFunctions.string(clinicalTrialJsonObject, "recruitment"))
                    .therapies(extractTherapyList(clinicalTrialJsonObject.getAsJsonArray("therapies")))
                    .build());
        }
        return clinicalTrials;
    }

    @NotNull
    private static List<TherapyInfo> extractTherapyList(@NotNull JsonArray jsonArray) {
        List<TherapyInfo> therapies = Lists.newArrayList();
        JsonDatamodelChecker therapyChecker = TherapyDataModelChecker.therapyChecker();

        for (JsonElement therapy : jsonArray) {
            JsonObject therapyJsonObject = therapy.getAsJsonObject();
            therapyChecker.check(therapyJsonObject);

            therapies.add(ImmutableTherapyInfo.builder()
                    .id(JsonFunctions.integer(therapyJsonObject, "id"))
                    .therapyName(JsonFunctions.string(therapyJsonObject, "therapyName"))
                    .synonyms(JsonFunctions.optionalStringList(therapyJsonObject, "synonyms"))
                    .build());
        }
        return therapies;
    }

    @NotNull
    private static List<DrugInfo> extractDrugs(@NotNull JsonArray jsonArray) {
        List<DrugInfo> drugs = Lists.newArrayList();
        JsonDatamodelChecker drugChecker = TherapyDataModelChecker.drugChecker();

        for (JsonElement drug : jsonArray) {
            JsonObject drugJsonObject = drug.getAsJsonObject();
            drugChecker.check(drugJsonObject);

            drugs.add(ImmutableDrugInfo.builder()
                    .id(JsonFunctions.integer(drugJsonObject, "id"))
                    .drugName(JsonFunctions.string(drugJsonObject, "drugName"))
                    .terms(JsonFunctions.stringList(drugJsonObject, "terms"))
                    .build());
        }
        return drugs;
    }

    @NotNull
    private static List<GlobalApprovalStatusInfo> extractGlobalApprovalStatuses(@NotNull JsonArray jsonArray) {
        List<GlobalApprovalStatusInfo> globalApprovalStatuses = Lists.newArrayList();
        JsonDatamodelChecker globalApprovalStatusChecker = TherapyDataModelChecker.globalApprovalStatusChecker();

        for (JsonElement globalApprovalStatus : jsonArray) {
            JsonObject globalApprovalStatusJsonObject = globalApprovalStatus.getAsJsonObject();
            globalApprovalStatusChecker.check(globalApprovalStatusJsonObject);

            globalApprovalStatuses.add(ImmutableGlobalApprovalStatusInfo.builder()
                    .id(JsonFunctions.integer(globalApprovalStatusJsonObject, "id"))
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