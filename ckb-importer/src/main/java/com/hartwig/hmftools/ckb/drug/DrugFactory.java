package com.hartwig.hmftools.ckb.drug;

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
import com.hartwig.hmftools.ckb.clinicaltrial.ClinicalTrialDataModelChecker;
import com.hartwig.hmftools.ckb.common.ImmutableIndicationInfo;
import com.hartwig.hmftools.ckb.common.ImmutableMolecularProfileInfo;
import com.hartwig.hmftools.ckb.common.ImmutableTherapyInfo;
import com.hartwig.hmftools.ckb.common.IndicationInfo;
import com.hartwig.hmftools.ckb.common.MolecularProfileInfo;
import com.hartwig.hmftools.ckb.common.TherapyInfo;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;
import com.hartwig.hmftools.common.utils.json.JsonFunctions;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class DrugFactory {

    private static final Logger LOGGER = LogManager.getLogger(DrugFactory.class);

    private DrugFactory() {
    }

    @NotNull
    public static List<Drug> readingDrugs(@NotNull String drugsDir) throws IOException {
        LOGGER.info("Start reading drug");

        List<Drug> drugs = Lists.newArrayList();
        File[] filesDrugs = new File(drugsDir).listFiles();

        if (filesDrugs != null) {
            LOGGER.info("The total files in the drug dir is {}", filesDrugs.length);

            for (File drug : filesDrugs) {
                JsonParser parser = new JsonParser();
                JsonReader reader = new JsonReader(new FileReader(drug));
                reader.setLenient(true);

                while (reader.peek() != JsonToken.END_DOCUMENT) {
                    JsonObject drugsEntryObject = parser.parse(reader).getAsJsonObject();
                    JsonDatamodelChecker drugChecker = DrugDataModelChecker.drugObjectChecker();
                    drugChecker.check(drugsEntryObject);

                    drugs.add(ImmutableDrug.builder()
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
                reader.close();
            }
        }
        LOGGER.info("Finished reading drug");

        return drugs;
    }

    @NotNull
    public static List<DrugDescription> extractDrugDescriptions(@NotNull JsonArray jsonArray) {
        List<DrugDescription> drugDescriptions = Lists.newArrayList();
        JsonDatamodelChecker drugDescriptionChecker = DrugDataModelChecker.drugDescriptionObjectChecker();
        for (JsonElement drugDescription : jsonArray) {
            JsonObject drugDescriptionObject = drugDescription.getAsJsonObject();
            drugDescriptionChecker.check(drugDescriptionObject);

            drugDescriptions.add(ImmutableDrugDescription.builder()
                    .description(JsonFunctions.string(drugDescriptionObject, "description"))
                    .references(extractDrugReferences(drugDescriptionObject.getAsJsonArray("references")))
                    .build());
        }
        return drugDescriptions;
    }

    @NotNull
    public static List<DrugReference> extractDrugReferences(@NotNull JsonArray jsonArray) {
        List<DrugReference> drugReferences = Lists.newArrayList();
        JsonDatamodelChecker drugReferenceChecker = DrugDataModelChecker.drugReferenceObjectChecker();

        for (JsonElement drugReference : jsonArray) {
            JsonObject drugsReferenceObject = drugReference.getAsJsonObject();
            drugReferenceChecker.check(drugsReferenceObject);

            drugReferences.add(ImmutableDrugReference.builder()
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
    public static List<DrugClass> extractDrugsClasses(@NotNull JsonArray jsonArray) {
        List<DrugClass> drugClasses = Lists.newArrayList();
        JsonDatamodelChecker drugClassChecker = DrugDataModelChecker.drugClassObjectChecker();

        for (JsonElement drugClass : jsonArray) {
            JsonObject drugClassObject = drugClass.getAsJsonObject();
            drugClassChecker.check(drugClassObject);

            drugClasses.add(ImmutableDrugClass.builder()
                    .id(JsonFunctions.string(drugClassObject, "id"))
                    .drugClass(JsonFunctions.string(drugClassObject, "drugClass"))
                    .build());
        }
        return drugClasses;
    }

    @NotNull
    public static List<DrugClinicalTrial> extractCliniclaTrials(@NotNull JsonArray jsonArray) {
        List<DrugClinicalTrial> clinicalTrials = Lists.newArrayList();
        JsonDatamodelChecker drugClinicalTrialChecker = DrugDataModelChecker.drugClinicalTrialObjectChecker();

        for (JsonElement clinicalTrial : jsonArray) {
            JsonObject clinicalTrialObject = clinicalTrial.getAsJsonObject();
            drugClinicalTrialChecker.check(clinicalTrialObject);

            clinicalTrials.add(ImmutableDrugClinicalTrial.builder()
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
    public static List<TherapyInfo> extractDrugTherapies(@NotNull JsonArray jsonArray) {
        List<TherapyInfo> drugTherapies = Lists.newArrayList();
        JsonDatamodelChecker drugClinicalTrialTherapyChecker = DrugDataModelChecker.drugClinicalTrialTherapyObjectChecker();

        for (JsonElement drugTherpy : jsonArray) {
            JsonObject drugTherapyObject = drugTherpy.getAsJsonObject();
            drugClinicalTrialTherapyChecker.check(drugTherapyObject);

            drugTherapies.add(ImmutableTherapyInfo.builder()
                    .id(JsonFunctions.string(drugTherapyObject, "id"))
                    .therapyName(JsonFunctions.string(drugTherapyObject, "therapyName"))
                    .synonyms(JsonFunctions.nullableString(drugTherapyObject, "synonyms"))
                    .build());
        }
        return drugTherapies;
    }

    @NotNull
    public static List<DrugEvidence> extractEvidence(@NotNull JsonArray jsonArray) {
        List<DrugEvidence> evidences = Lists.newArrayList();
        JsonDatamodelChecker drugEvidenceChecker = DrugDataModelChecker.drugEvidenceObjectChecker();

        for (JsonElement evidence : jsonArray) {
            JsonObject evidenceObject = evidence.getAsJsonObject();
            drugEvidenceChecker.check(evidenceObject);

            evidences.add(ImmutableDrugEvidence.builder()
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
    public static MolecularProfileInfo extractMolecularProfile(@NotNull JsonObject jsonObject) {
        JsonObject molecularProfileObject = jsonObject.getAsJsonObject();
        JsonDatamodelChecker drugEvidenceMolecularProfileChecker = DrugDataModelChecker.drugEvidenceMolecularProfileObjectChecker();
        drugEvidenceMolecularProfileChecker.check(molecularProfileObject);

        return ImmutableMolecularProfileInfo.builder()
                .id(JsonFunctions.string(molecularProfileObject, "id"))
                .profileName(JsonFunctions.string(molecularProfileObject, "profileName"))
                .profileTreatmentApproache(null)
                .build();
    }

    @NotNull
    public static TherapyInfo extractTherapy(@NotNull JsonObject jsonObject) {
        JsonObject therapyObject = jsonObject.getAsJsonObject();
        JsonDatamodelChecker drugEvidenceTherapyChecker = DrugDataModelChecker.drugEvidenceTherapyObjectChecker();
        drugEvidenceTherapyChecker.check(therapyObject);

        return ImmutableTherapyInfo.builder()
                .id(JsonFunctions.string(therapyObject, "id"))
                .therapyName(JsonFunctions.string(therapyObject, "therapyName"))
                .synonyms(JsonFunctions.nullableString(therapyObject, "synonyms"))
                .build();
    }

    @NotNull
    public static IndicationInfo extractIndications(@NotNull JsonObject jsonObject) {
        JsonObject indicationObject = jsonObject.getAsJsonObject();
        JsonDatamodelChecker drugEvidenceIndicationChecker = DrugDataModelChecker.drugEvidenceIndicationObjectChecker();
        drugEvidenceIndicationChecker.check(indicationObject);

        return ImmutableIndicationInfo.builder()
                .id(JsonFunctions.string(indicationObject, "id"))
                .name(JsonFunctions.string(indicationObject, "name"))
                .source(JsonFunctions.string(indicationObject, "source"))
                .build();
    }

    @NotNull
    public static List<DrugEvidenceReference> extractEvidenceReferences(@NotNull JsonArray jsonArray) {
        List<DrugEvidenceReference> references = Lists.newArrayList();
        JsonDatamodelChecker drugEvidenceReferenceChecker = DrugDataModelChecker.drugEvidenceReferenceObjectChecker();

        for (JsonElement reference : jsonArray) {
            JsonObject referenceObject = reference.getAsJsonObject();
            drugEvidenceReferenceChecker.check(referenceObject);

            references.add(ImmutableDrugEvidenceReference.builder()
                    .id(JsonFunctions.string(referenceObject, "id"))
                    .pubMedId(JsonFunctions.nullableString(referenceObject, "pubMedId"))
                    .title(JsonFunctions.nullableString(referenceObject, "title"))
                    .url(JsonFunctions.nullableString(referenceObject, "url"))
                    .build());
        }
        return references;
    }

    @NotNull
    public static List<TherapyInfo> extractTherapies(@NotNull JsonArray jsonArray) {
        List<TherapyInfo> drugsTherapies = Lists.newArrayList();
        JsonDatamodelChecker drugTherapyChecker = DrugDataModelChecker.drugTherapyObjectChecker();

        for (JsonElement drugTherapy : jsonArray) {
            JsonObject drugTherapyObject = drugTherapy.getAsJsonObject();
            drugTherapyChecker.check(drugTherapyObject);

            drugsTherapies.add(ImmutableTherapyInfo.builder()
                    .id(JsonFunctions.string(drugTherapyObject, "id"))
                    .therapyName(JsonFunctions.string(drugTherapyObject, "therapyName"))
                    .synonyms(JsonFunctions.nullableString(drugTherapyObject, "synonyms"))
                    .build());
        }
        return drugsTherapies;
    }

    @NotNull
    public static List<DrugGlobalApprovalStatus> extractGlobalApprovaStatus(@NotNull JsonArray jsonArray) {
        List<DrugGlobalApprovalStatus> globalApproavalStatuses = Lists.newArrayList();
        JsonDatamodelChecker drugGlobalApprovalStatusChecker = DrugDataModelChecker.drugGlobalApprovalStatusObjectChecker();

        for (JsonElement globalApproavalStatus : jsonArray) {
            JsonObject globalTherapyObject = globalApproavalStatus.getAsJsonObject();
            drugGlobalApprovalStatusChecker.check(globalTherapyObject);

            globalApproavalStatuses.add(ImmutableDrugGlobalApprovalStatus.builder()
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
    public static TherapyInfo extractTherapiesGlobalApprovalStatus(@NotNull JsonObject jsonObject) {
        JsonObject therapiesGlobalApprovalStatusObject = jsonObject.getAsJsonObject();
        JsonDatamodelChecker drugGlobalApprovalStatusTherapyChecker = DrugDataModelChecker.drugGlobalApprovalStatusTherapyObjectChecker();
        drugGlobalApprovalStatusTherapyChecker.check(therapiesGlobalApprovalStatusObject);

        return ImmutableTherapyInfo.builder()
                .id(JsonFunctions.string(therapiesGlobalApprovalStatusObject, "id"))
                .therapyName(JsonFunctions.string(therapiesGlobalApprovalStatusObject, "therapyName"))
                .synonyms(JsonFunctions.string(therapiesGlobalApprovalStatusObject, "synonyms"))
                .build();
    }

    @NotNull
    public static IndicationInfo extractIndicationsGlobalApprovalStatus(@NotNull JsonObject jsonObject) {
        JsonObject indicationsGlobalApprovalStatusObject = jsonObject.getAsJsonObject();
        JsonDatamodelChecker drugGlobalApprovalStatusIndicationChecker = DrugDataModelChecker.drugGlobalApprovalStatusIndicationObjectChecker();
        drugGlobalApprovalStatusIndicationChecker.check(indicationsGlobalApprovalStatusObject);

        return ImmutableIndicationInfo.builder()
                .id(JsonFunctions.string(indicationsGlobalApprovalStatusObject, "id"))
                .name(JsonFunctions.string(indicationsGlobalApprovalStatusObject, "name"))
                .source(JsonFunctions.string(indicationsGlobalApprovalStatusObject, "source"))
                .build();
    }

    @NotNull
    public static MolecularProfileInfo extractMolecularProfileGlobalApprovalStatus(@NotNull JsonObject jsonObject) {
        JsonObject molecularProfileGlobalApprovalStatus = jsonObject.getAsJsonObject();
        JsonDatamodelChecker drugGlobalApprovalStatusMolecularProfileChecker = DrugDataModelChecker.drugGlobalApprovalStatusMolecularProfileObjectChecker();
        drugGlobalApprovalStatusMolecularProfileChecker.check(molecularProfileGlobalApprovalStatus);

        return ImmutableMolecularProfileInfo.builder()
                .id(JsonFunctions.string(molecularProfileGlobalApprovalStatus, "id"))
                .profileName(JsonFunctions.string(molecularProfileGlobalApprovalStatus, "profileName"))
                .build();
    }
}
