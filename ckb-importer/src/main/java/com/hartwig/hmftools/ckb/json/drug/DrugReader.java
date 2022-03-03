package com.hartwig.hmftools.ckb.json.drug;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.ckb.json.CkbJsonDirectoryReader;
import com.hartwig.hmftools.ckb.json.common.ClinicalTrialInfo;
import com.hartwig.hmftools.ckb.json.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.json.common.DrugClassInfo;
import com.hartwig.hmftools.ckb.json.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.json.common.GlobalApprovalStatusInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableClinicalTrialInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableDescriptionInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableDrugClassInfo;
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

public class DrugReader extends CkbJsonDirectoryReader<JsonDrug> {
    private static final Logger LOGGER = LogManager.getLogger(DrugReader.class);

    public DrugReader(@Nullable final Integer maxFilesToRead) {
        super(maxFilesToRead);
    }

    @NotNull
    @Override
    protected JsonDrug read(@NotNull final JsonObject object) {
        JsonDatamodelChecker drugChecker = DrugDataModelChecker.drugObjectChecker();
        drugChecker.check(object);

        return ImmutableJsonDrug.builder()
                .id(JsonFunctions.integer(object, "id"))
                .drugName(JsonFunctions.string(object, "drugName"))
                .terms(JsonFunctions.stringList(object, "terms"))
                .synonyms(JsonFunctions.stringList(object, "synonyms"))
                .tradeName(JsonFunctions.nullableString(object, "tradeName"))
                .descriptions(extractDescriptions(object.getAsJsonArray("drugDescriptions")))
                .drugClasses(extractDrugsClasses(object.getAsJsonArray("drugClasses")))
                .casRegistryNum(JsonFunctions.nullableString(object, "casRegistryNum"))
                .ncitId(JsonFunctions.nullableString(object, "ncitId"))
                .createDate(DateConverter.toDate(JsonFunctions.string(object, "createDate")))
                .clinicalTrials(extractClinicalTrials(object.getAsJsonArray("clinicalTrials")))
                .evidence(extractEvidence(object.getAsJsonArray("evidence")))
                .therapies(extractTherapies(object.getAsJsonArray("therapies")))
                .globalApprovalStatus(extractGlobalApprovalStatus(object.getAsJsonArray("globalApprovalStatus")))
                .build();
    }

    @NotNull
    private static List<DescriptionInfo> extractDescriptions(@NotNull JsonArray jsonArray) {
        List<DescriptionInfo> descriptions = Lists.newArrayList();
        JsonDatamodelChecker descriptionChecker = DrugDataModelChecker.descriptionObjectChecker();
        for (JsonElement description : jsonArray) {
            JsonObject descriptionObject = description.getAsJsonObject();
            descriptionChecker.check(descriptionObject);

            descriptions.add(ImmutableDescriptionInfo.builder()
                    .description(JsonFunctions.string(descriptionObject, "description"))
                    .references(extractDrugReferences(descriptionObject.getAsJsonArray("references")))
                    .build());
        }
        return descriptions;
    }

    @NotNull
    private static List<ReferenceInfo> extractDrugReferences(@NotNull JsonArray jsonArray) {
        List<ReferenceInfo> references = Lists.newArrayList();
        JsonDatamodelChecker referenceChecker = DrugDataModelChecker.referenceObjectChecker();

        for (JsonElement reference : jsonArray) {
            JsonObject referenceObject = reference.getAsJsonObject();
            referenceChecker.check(referenceObject);

            references.add(ImmutableReferenceInfo.builder()
                    .id(JsonFunctions.integer(referenceObject, "id"))
                    .pubMedId(JsonFunctions.nullableString(referenceObject, "pubMedId"))
                    .title(JsonFunctions.nullableString(referenceObject, "title"))
                    .url(JsonFunctions.nullableString(referenceObject, "url"))
                    .authors(JsonFunctions.nullableString(referenceObject, "authors"))
                    .journal(JsonFunctions.nullableString(referenceObject, "journal"))
                    .volume(JsonFunctions.nullableString(referenceObject, "volume"))
                    .issue(JsonFunctions.nullableString(referenceObject, "issue"))
                    .date(JsonFunctions.nullableString(referenceObject, "date"))
                    .abstractText(JsonFunctions.nullableString(referenceObject, "abstractText"))
                    .year(JsonFunctions.nullableString(referenceObject, "year"))
                    .build());
        }
        return references;
    }

    @NotNull
    private static List<DrugClassInfo> extractDrugsClasses(@NotNull JsonArray jsonArray) {
        List<DrugClassInfo> drugClasses = Lists.newArrayList();
        JsonDatamodelChecker drugClassChecker = DrugDataModelChecker.drugClassObjectChecker();

        for (JsonElement drugClass : jsonArray) {
            JsonObject drugClassObject = drugClass.getAsJsonObject();
            drugClassChecker.check(drugClassObject);

            drugClasses.add(ImmutableDrugClassInfo.builder()
                    .id(JsonFunctions.integer(drugClassObject, "id"))
                    .drugClass(JsonFunctions.string(drugClassObject, "drugClass"))
                    .build());
        }
        return drugClasses;
    }

    @NotNull
    private static List<ClinicalTrialInfo> extractClinicalTrials(@NotNull JsonArray jsonArray) {
        List<ClinicalTrialInfo> clinicalTrials = Lists.newArrayList();
        JsonDatamodelChecker clinicalTrialChecker = DrugDataModelChecker.clinicalTrialObjectChecker();

        for (JsonElement clinicalTrial : jsonArray) {
            JsonObject clinicalTrialObject = clinicalTrial.getAsJsonObject();
            clinicalTrialChecker.check(clinicalTrialObject);

            if (JsonFunctions.nullableString(clinicalTrialObject, "phase") == null) {
                LOGGER.warn("phase of study '{}' is nullable from DrugReader", JsonFunctions.string(clinicalTrialObject, "nctId"));
            }

            clinicalTrials.add(ImmutableClinicalTrialInfo.builder()
                    .nctId(JsonFunctions.string(clinicalTrialObject, "nctId"))
                    .title(JsonFunctions.string(clinicalTrialObject, "title"))
                    .phase(JsonFunctions.nullableString(clinicalTrialObject, "phase"))
                    .recruitment(JsonFunctions.string(clinicalTrialObject, "recruitment"))
                    .therapies(extractClinicalTrialTherapies(clinicalTrialObject.getAsJsonArray("therapies")))
                    .build());
        }
        return clinicalTrials;
    }

    @NotNull
    private static List<TherapyInfo> extractClinicalTrialTherapies(@NotNull JsonArray jsonArray) {
        List<TherapyInfo> therapies = Lists.newArrayList();
        JsonDatamodelChecker clinicalTrialTherapyChecker = DrugDataModelChecker.clinicalTrialTherapyObjectChecker();

        for (JsonElement therapy : jsonArray) {
            JsonObject therapyObject = therapy.getAsJsonObject();
            clinicalTrialTherapyChecker.check(therapyObject);

            therapies.add(ImmutableTherapyInfo.builder()
                    .id(JsonFunctions.integer(therapyObject, "id"))
                    .therapyName(JsonFunctions.string(therapyObject, "therapyName"))
                    .synonyms(JsonFunctions.optionalStringList(therapyObject, "synonyms"))
                    .build());
        }
        return therapies;
    }

    @NotNull
    private static List<EvidenceInfo> extractEvidence(@NotNull JsonArray jsonArray) {
        List<EvidenceInfo> evidences = Lists.newArrayList();
        JsonDatamodelChecker evidenceChecker = DrugDataModelChecker.evidenceObjectChecker();

        for (JsonElement evidence : jsonArray) {
            JsonObject evidenceObject = evidence.getAsJsonObject();
            evidenceChecker.check(evidenceObject);

            evidences.add(ImmutableEvidenceInfo.builder()
                    .id(JsonFunctions.integer(evidenceObject, "id"))
                    .approvalStatus(JsonFunctions.string(evidenceObject, "approvalStatus"))
                    .evidenceType(JsonFunctions.string(evidenceObject, "evidenceType"))
                    .efficacyEvidence(JsonFunctions.string(evidenceObject, "efficacyEvidence"))
                    .molecularProfile(extractMolecularProfile(evidenceObject.getAsJsonObject("molecularProfile")))
                    .therapy(extractTherapy(evidenceObject.getAsJsonObject("therapy")))
                    .indication(extractIndication(evidenceObject.getAsJsonObject("indication")))
                    .responseType(JsonFunctions.string(evidenceObject, "responseType"))
                    .references(extractEvidenceReferences(evidenceObject.getAsJsonArray("references")))
                    .ampCapAscoEvidenceLevel(JsonFunctions.string(evidenceObject, "ampCapAscoEvidenceLevel"))
                    .ampCapAscoInferredTier(JsonFunctions.string(evidenceObject, "ampCapAscoInferredTier"))
                    .build());
        }
        return evidences;
    }

    @NotNull
    private static MolecularProfileInfo extractMolecularProfile(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker evidenceMolecularProfileChecker = DrugDataModelChecker.evidenceMolecularProfileObjectChecker();
        evidenceMolecularProfileChecker.check(jsonObject);

        return ImmutableMolecularProfileInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .profileName(JsonFunctions.string(jsonObject, "profileName"))
                .build();
    }

    @NotNull
    private static TherapyInfo extractTherapy(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker evidenceTherapyChecker = DrugDataModelChecker.evidenceTherapyObjectChecker();
        evidenceTherapyChecker.check(jsonObject);

        return ImmutableTherapyInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .therapyName(JsonFunctions.string(jsonObject, "therapyName"))
                .synonyms(JsonFunctions.optionalStringList(jsonObject, "synonyms"))
                .build();
    }

    @NotNull
    private static IndicationInfo extractIndication(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker evidenceIndicationChecker = DrugDataModelChecker.evidenceIndicationObjectChecker();
        evidenceIndicationChecker.check(jsonObject);

        return ImmutableIndicationInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .name(JsonFunctions.string(jsonObject, "name"))
                .source(JsonFunctions.string(jsonObject, "source"))
                .build();
    }

    @NotNull
    private static List<ReferenceInfo> extractEvidenceReferences(@NotNull JsonArray jsonArray) {
        List<ReferenceInfo> references = Lists.newArrayList();
        JsonDatamodelChecker evidenceReferenceChecker = DrugDataModelChecker.evidenceReferenceObjectChecker();

        for (JsonElement reference : jsonArray) {
            JsonObject referenceObject = reference.getAsJsonObject();
            evidenceReferenceChecker.check(referenceObject);

            references.add(ImmutableReferenceInfo.builder()
                    .id(JsonFunctions.integer(referenceObject, "id"))
                    .pubMedId(JsonFunctions.nullableString(referenceObject, "pubMedId"))
                    .title(JsonFunctions.nullableString(referenceObject, "title"))
                    .url(JsonFunctions.nullableString(referenceObject, "url"))
                    .build());
        }
        return references;
    }

    @NotNull
    private static List<TherapyInfo> extractTherapies(@NotNull JsonArray jsonArray) {
        List<TherapyInfo> therapies = Lists.newArrayList();
        JsonDatamodelChecker therapyChecker = DrugDataModelChecker.therapyObjectChecker();

        for (JsonElement therapy : jsonArray) {
            JsonObject therapyObject = therapy.getAsJsonObject();
            therapyChecker.check(therapyObject);

            therapies.add(ImmutableTherapyInfo.builder()
                    .id(JsonFunctions.integer(therapyObject, "id"))
                    .therapyName(JsonFunctions.string(therapyObject, "therapyName"))
                    .synonyms(JsonFunctions.optionalStringList(therapyObject, "synonyms"))
                    .build());
        }
        return therapies;
    }

    @NotNull
    private static List<GlobalApprovalStatusInfo> extractGlobalApprovalStatus(@NotNull JsonArray jsonArray) {
        List<GlobalApprovalStatusInfo> globalApprovalStatuses = Lists.newArrayList();
        JsonDatamodelChecker drugGlobalApprovalStatusChecker = DrugDataModelChecker.globalApprovalStatusObjectChecker();

        for (JsonElement globalApprovalStatus : jsonArray) {
            JsonObject globalStatusObject = globalApprovalStatus.getAsJsonObject();
            drugGlobalApprovalStatusChecker.check(globalStatusObject);

            globalApprovalStatuses.add(ImmutableGlobalApprovalStatusInfo.builder()
                    .id(JsonFunctions.integer(globalStatusObject, "id"))
                    .therapy(extractTherapyGlobalApprovalStatus(globalStatusObject.getAsJsonObject("therapy")))
                    .indication(extractIndicationGlobalApprovalStatus(globalStatusObject.getAsJsonObject("indication")))
                    .molecularProfile(extractMolecularProfileGlobalApprovalStatus(globalStatusObject.getAsJsonObject("molecularProfile")))
                    .approvalAuthority(JsonFunctions.string(globalStatusObject, "approvalAuthority"))
                    .approvalStatus(JsonFunctions.string(globalStatusObject, "approvalStatus"))
                    .build());
        }
        return globalApprovalStatuses;
    }

    @NotNull
    private static TherapyInfo extractTherapyGlobalApprovalStatus(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker globalApprovalStatusTherapyChecker = DrugDataModelChecker.globalApprovalStatusTherapyObjectChecker();
        globalApprovalStatusTherapyChecker.check(jsonObject);

        return ImmutableTherapyInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .therapyName(JsonFunctions.string(jsonObject, "therapyName"))
                .synonyms(JsonFunctions.optionalStringList(jsonObject, "synonyms"))
                .build();
    }

    @NotNull
    private static IndicationInfo extractIndicationGlobalApprovalStatus(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker globalApprovalStatusIndicationChecker = DrugDataModelChecker.globalApprovalStatusIndicationObjectChecker();
        globalApprovalStatusIndicationChecker.check(jsonObject);

        return ImmutableIndicationInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .name(JsonFunctions.string(jsonObject, "name"))
                .source(JsonFunctions.string(jsonObject, "source"))
                .build();
    }

    @NotNull
    private static MolecularProfileInfo extractMolecularProfileGlobalApprovalStatus(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker globalApprovalStatusMolecularProfileChecker =
                DrugDataModelChecker.globalApprovalStatusMolecularProfileObjectChecker();
        globalApprovalStatusMolecularProfileChecker.check(jsonObject);

        return ImmutableMolecularProfileInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .profileName(JsonFunctions.string(jsonObject, "profileName"))
                .build();
    }
}
