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

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class DrugReader extends CkbJsonDirectoryReader<Drug> {

    public DrugReader(@Nullable final Integer maxFilesToRead) {
        super(maxFilesToRead);
    }

    @NotNull
    @Override
    protected Drug read(@NotNull final JsonObject object) {
        JsonDatamodelChecker drugChecker = DrugDataModelChecker.drugObjectChecker();
        drugChecker.check(object);

        return ImmutableDrug.builder()
                .id(JsonFunctions.integer(object, "id"))
                .drugName(JsonFunctions.string(object, "drugName"))
                .term(JsonFunctions.stringList(object, "terms"))
                .synonym(JsonFunctions.stringList(object, "synonyms"))
                .tradeName(JsonFunctions.nullableString(object, "tradeName"))
                .description(extractDrugDescriptions(object.getAsJsonArray("drugDescriptions")))
                .drugClass(extractDrugsClasses(object.getAsJsonArray("drugClasses")))
                .casRegistryNum(JsonFunctions.nullableString(object, "casRegistryNum"))
                .nctId(JsonFunctions.nullableString(object, "ncitId"))
                .createDate(DateConverter.toDate(JsonFunctions.string(object, "createDate")))
                .clinicalTrial(extractCliniclaTrials(object.getAsJsonArray("clinicalTrials")))
                .evidence(extractEvidence(object.getAsJsonArray("evidence")))
                .therapy(extractTherapies(object.getAsJsonArray("therapies")))
                .globalApprovalStatus(object.has("globalApprovaStatus") ? extractGlobalApprovaStatus(object.getAsJsonArray(
                        "globalApprovaStatus")) : null)
                .build();
    }

    @NotNull
    private static List<DescriptionInfo> extractDrugDescriptions(@NotNull JsonArray jsonArray) {
        List<DescriptionInfo> drugDescriptions = Lists.newArrayList();
        JsonDatamodelChecker drugDescriptionChecker = DrugDataModelChecker.drugDescriptionObjectChecker();
        for (JsonElement drugDescription : jsonArray) {
            JsonObject drugDescriptionObject = drugDescription.getAsJsonObject();
            drugDescriptionChecker.check(drugDescriptionObject);

            drugDescriptions.add(ImmutableDescriptionInfo.builder()
                    .description(JsonFunctions.string(drugDescriptionObject, "description"))
                    .reference(extractDrugReferences(drugDescriptionObject.getAsJsonArray("references")))
                    .build());
        }
        return drugDescriptions;
    }

    @NotNull
    private static List<ReferenceInfo> extractDrugReferences(@NotNull JsonArray jsonArray) {
        List<ReferenceInfo> drugReferences = Lists.newArrayList();
        JsonDatamodelChecker drugReferenceChecker = DrugDataModelChecker.drugReferenceObjectChecker();

        for (JsonElement drugReference : jsonArray) {
            JsonObject drugsReferenceObject = drugReference.getAsJsonObject();
            drugReferenceChecker.check(drugsReferenceObject);

            drugReferences.add(ImmutableReferenceInfo.builder()
                    .id(JsonFunctions.integer(drugsReferenceObject, "id"))
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
    private static List<ClinicalTrialInfo> extractCliniclaTrials(@NotNull JsonArray jsonArray) {
        List<ClinicalTrialInfo> clinicalTrials = Lists.newArrayList();
        JsonDatamodelChecker drugClinicalTrialChecker = DrugDataModelChecker.drugClinicalTrialObjectChecker();

        for (JsonElement clinicalTrial : jsonArray) {
            JsonObject clinicalTrialObject = clinicalTrial.getAsJsonObject();
            drugClinicalTrialChecker.check(clinicalTrialObject);

            clinicalTrials.add(ImmutableClinicalTrialInfo.builder()
                    .nctId(JsonFunctions.string(clinicalTrialObject, "nctId"))
                    .title(JsonFunctions.string(clinicalTrialObject, "title"))
                    .phase(JsonFunctions.string(clinicalTrialObject, "phase"))
                    .recruitment(JsonFunctions.string(clinicalTrialObject, "recruitment"))
                    .therapy(extractDrugTherapies(clinicalTrialObject.getAsJsonArray("therapies")))
                    .build());
        }
        return clinicalTrials;
    }

    @NotNull
    private static List<TherapyInfo> extractDrugTherapies(@NotNull JsonArray jsonArray) {
        List<TherapyInfo> drugTherapies = Lists.newArrayList();
        JsonDatamodelChecker drugClinicalTrialTherapyChecker = DrugDataModelChecker.drugClinicalTrialTherapyObjectChecker();

        for (JsonElement drugTherpy : jsonArray) {
            JsonObject drugTherapyObject = drugTherpy.getAsJsonObject();
            drugClinicalTrialTherapyChecker.check(drugTherapyObject);

            drugTherapies.add(ImmutableTherapyInfo.builder()
                    .id(JsonFunctions.integer(drugTherapyObject, "id"))
                    .therapyName(JsonFunctions.string(drugTherapyObject, "therapyName"))
                    .synonyms(JsonFunctions.nullableString(drugTherapyObject, "synonyms"))
                    .build());
        }
        return drugTherapies;
    }

    @NotNull
    private static List<EvidenceInfo> extractEvidence(@NotNull JsonArray jsonArray) {
        List<EvidenceInfo> evidences = Lists.newArrayList();
        JsonDatamodelChecker drugEvidenceChecker = DrugDataModelChecker.drugEvidenceObjectChecker();

        for (JsonElement evidence : jsonArray) {
            JsonObject evidenceObject = evidence.getAsJsonObject();
            drugEvidenceChecker.check(evidenceObject);

            evidences.add(ImmutableEvidenceInfo.builder()
                    .id(JsonFunctions.integer(evidenceObject, "id"))
                    .approvalStatus(JsonFunctions.string(evidenceObject, "approvalStatus"))
                    .evidenceType(JsonFunctions.string(evidenceObject, "evidenceType"))
                    .efficacyEvidence(JsonFunctions.string(evidenceObject, "efficacyEvidence"))
                    .molecularProfile(extractMolecularProfile(evidenceObject.getAsJsonObject("molecularProfile")))
                    .therapy(extractTherapy(evidenceObject.getAsJsonObject("therapy")))
                    .indication(evidenceObject.has("indications")
                            ? extractIndications(evidenceObject.getAsJsonObject("indications"))
                            : null)
                    .responseType(JsonFunctions.string(evidenceObject, "responseType"))
                    .reference(extractEvidenceReferences(evidenceObject.getAsJsonArray("references")))
                    .ampCapAscoEvidenceLevel(JsonFunctions.string(evidenceObject, "ampCapAscoEvidenceLevel"))
                    .ampCapAscoInferredTier(JsonFunctions.string(evidenceObject, "ampCapAscoInferredTier"))
                    .build());
        }
        return evidences;
    }

    @NotNull
    private static MolecularProfileInfo extractMolecularProfile(@NotNull JsonObject jsonObject) {
        JsonObject molecularProfileObject = jsonObject.getAsJsonObject();
        JsonDatamodelChecker drugEvidenceMolecularProfileChecker = DrugDataModelChecker.drugEvidenceMolecularProfileObjectChecker();
        drugEvidenceMolecularProfileChecker.check(molecularProfileObject);

        return ImmutableMolecularProfileInfo.builder()
                .id(JsonFunctions.integer(molecularProfileObject, "id"))
                .profileName(JsonFunctions.string(molecularProfileObject, "profileName"))
                .build();
    }

    @NotNull
    private static TherapyInfo extractTherapy(@NotNull JsonObject jsonObject) {
        JsonObject therapyObject = jsonObject.getAsJsonObject();
        JsonDatamodelChecker drugEvidenceTherapyChecker = DrugDataModelChecker.drugEvidenceTherapyObjectChecker();
        drugEvidenceTherapyChecker.check(therapyObject);

        return ImmutableTherapyInfo.builder()
                .id(JsonFunctions.integer(therapyObject, "id"))
                .therapyName(JsonFunctions.string(therapyObject, "therapyName"))
                .synonyms(JsonFunctions.nullableString(therapyObject, "synonyms"))
                .build();
    }

    @NotNull
    private static IndicationInfo extractIndications(@NotNull JsonObject jsonObject) {
        JsonObject indicationObject = jsonObject.getAsJsonObject();
        JsonDatamodelChecker drugEvidenceIndicationChecker = DrugDataModelChecker.drugEvidenceIndicationObjectChecker();
        drugEvidenceIndicationChecker.check(indicationObject);

        return ImmutableIndicationInfo.builder()
                .id(JsonFunctions.integer(indicationObject, "id"))
                .name(JsonFunctions.string(indicationObject, "name"))
                .source(JsonFunctions.string(indicationObject, "source"))
                .build();
    }

    @NotNull
    private static List<ReferenceInfo> extractEvidenceReferences(@NotNull JsonArray jsonArray) {
        List<ReferenceInfo> references = Lists.newArrayList();
        JsonDatamodelChecker drugEvidenceReferenceChecker = DrugDataModelChecker.drugEvidenceReferenceObjectChecker();

        for (JsonElement reference : jsonArray) {
            JsonObject referenceObject = reference.getAsJsonObject();
            drugEvidenceReferenceChecker.check(referenceObject);

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
        List<TherapyInfo> drugsTherapies = Lists.newArrayList();
        JsonDatamodelChecker drugTherapyChecker = DrugDataModelChecker.drugTherapyObjectChecker();

        for (JsonElement drugTherapy : jsonArray) {
            JsonObject drugTherapyObject = drugTherapy.getAsJsonObject();
            drugTherapyChecker.check(drugTherapyObject);

            drugsTherapies.add(ImmutableTherapyInfo.builder()
                    .id(JsonFunctions.integer(drugTherapyObject, "id"))
                    .therapyName(JsonFunctions.string(drugTherapyObject, "therapyName"))
                    .synonyms(JsonFunctions.nullableString(drugTherapyObject, "synonyms"))
                    .build());
        }
        return drugsTherapies;
    }

    @NotNull
    private static List<GlobalApprovalStatusInfo> extractGlobalApprovaStatus(@NotNull JsonArray jsonArray) {
        List<GlobalApprovalStatusInfo> globalApproavalStatuses = Lists.newArrayList();
        JsonDatamodelChecker drugGlobalApprovalStatusChecker = DrugDataModelChecker.drugGlobalApprovalStatusObjectChecker();

        for (JsonElement globalApproavalStatus : jsonArray) {
            JsonObject globalTherapyObject = globalApproavalStatus.getAsJsonObject();
            drugGlobalApprovalStatusChecker.check(globalTherapyObject);

            globalApproavalStatuses.add(ImmutableGlobalApprovalStatusInfo.builder()
                    .id(JsonFunctions.integer(globalTherapyObject, "id"))
                    .therapy(extractTherapiesGlobalApprovalStatus(globalTherapyObject.getAsJsonObject("therapy")))
                    .indication(extractIndicationsGlobalApprovalStatus(globalTherapyObject.getAsJsonObject("indications")))
                    .molecularProfile(extractMolecularProfileGlobalApprovalStatus(globalTherapyObject.getAsJsonObject("molecularProfile")))
                    .approvalAuthority(JsonFunctions.string(globalTherapyObject, "approvalAuthority"))
                    .approvalStatus(JsonFunctions.string(globalTherapyObject, "approvalStatus"))
                    .build());
        }
        return globalApproavalStatuses;
    }

    @NotNull
    private static TherapyInfo extractTherapiesGlobalApprovalStatus(@NotNull JsonObject jsonObject) {
        JsonObject therapiesGlobalApprovalStatusObject = jsonObject.getAsJsonObject();
        JsonDatamodelChecker drugGlobalApprovalStatusTherapyChecker = DrugDataModelChecker.drugGlobalApprovalStatusTherapyObjectChecker();
        drugGlobalApprovalStatusTherapyChecker.check(therapiesGlobalApprovalStatusObject);

        return ImmutableTherapyInfo.builder()
                .id(JsonFunctions.integer(therapiesGlobalApprovalStatusObject, "id"))
                .therapyName(JsonFunctions.string(therapiesGlobalApprovalStatusObject, "therapyName"))
                .synonyms(JsonFunctions.string(therapiesGlobalApprovalStatusObject, "synonyms"))
                .build();
    }

    @NotNull
    private static IndicationInfo extractIndicationsGlobalApprovalStatus(@NotNull JsonObject jsonObject) {
        JsonObject indicationsGlobalApprovalStatusObject = jsonObject.getAsJsonObject();
        JsonDatamodelChecker drugGlobalApprovalStatusIndicationChecker =
                DrugDataModelChecker.drugGlobalApprovalStatusIndicationObjectChecker();
        drugGlobalApprovalStatusIndicationChecker.check(indicationsGlobalApprovalStatusObject);

        return ImmutableIndicationInfo.builder()
                .id(JsonFunctions.integer(indicationsGlobalApprovalStatusObject, "id"))
                .name(JsonFunctions.string(indicationsGlobalApprovalStatusObject, "name"))
                .source(JsonFunctions.string(indicationsGlobalApprovalStatusObject, "source"))
                .build();
    }

    @NotNull
    private static MolecularProfileInfo extractMolecularProfileGlobalApprovalStatus(@NotNull JsonObject jsonObject) {
        JsonObject molecularProfileGlobalApprovalStatus = jsonObject.getAsJsonObject();
        JsonDatamodelChecker drugGlobalApprovalStatusMolecularProfileChecker =
                DrugDataModelChecker.drugGlobalApprovalStatusMolecularProfileObjectChecker();
        drugGlobalApprovalStatusMolecularProfileChecker.check(molecularProfileGlobalApprovalStatus);

        return ImmutableMolecularProfileInfo.builder()
                .id(JsonFunctions.integer(molecularProfileGlobalApprovalStatus, "id"))
                .profileName(JsonFunctions.string(molecularProfileGlobalApprovalStatus, "profileName"))
                .build();
    }
}
