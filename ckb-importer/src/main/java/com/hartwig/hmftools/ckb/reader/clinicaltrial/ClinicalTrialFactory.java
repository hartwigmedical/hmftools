package com.hartwig.hmftools.ckb.reader.clinicaltrial;

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

import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrial;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrialContact;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrialLocation;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrialVariantRequirementDetail;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ImmutableClinicalTrial;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ImmutableClinicalTrialContact;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ImmutableClinicalTrialLocation;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ImmutableClinicalTrialVariantRequirementDetail;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableIndicationInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableMolecularProfileInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableTherapyInfo;
import com.hartwig.hmftools.ckb.datamodel.common.IndicationInfo;
import com.hartwig.hmftools.ckb.datamodel.common.MolecularProfileInfo;
import com.hartwig.hmftools.ckb.datamodel.common.TherapyInfo;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;
import com.hartwig.hmftools.common.utils.json.JsonFunctions;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ClinicalTrialFactory {

    private static final Logger LOGGER = LogManager.getLogger(ClinicalTrialFactory.class);

    private ClinicalTrialFactory() {

    }

    @NotNull
    public static List<ClinicalTrial> readingClinicalTrial(@NotNull String clinicalTrialDir) throws IOException {
        LOGGER.info("Start reading clinical trials");

        List<ClinicalTrial> clinicalTrials = Lists.newArrayList();
        File[] filesClinicalTrials = new File(clinicalTrialDir).listFiles();

        if (filesClinicalTrials != null) {
            LOGGER.info("The total files in the clinical trial dir is {}", filesClinicalTrials.length);

            for (File clinicalTrial : filesClinicalTrials) {
                JsonParser parser = new JsonParser();
                JsonReader reader = new JsonReader(new FileReader(clinicalTrial));
                reader.setLenient(true);

                while (reader.peek() != JsonToken.END_DOCUMENT) {
                    JsonObject clinicalTrialsEntryObject = parser.parse(reader).getAsJsonObject();
                    JsonDatamodelChecker clinicalTrailChecker = ClinicalTrialDataModelChecker.clinicalTrialObjectChecker();
                    clinicalTrailChecker.check(clinicalTrialsEntryObject);

                    clinicalTrials.add(ImmutableClinicalTrial.builder()
                            .nctId(JsonFunctions.string(clinicalTrialsEntryObject, "nctId"))
                            .title(JsonFunctions.string(clinicalTrialsEntryObject, "title"))
                            .phase(JsonFunctions.string(clinicalTrialsEntryObject, "phase"))
                            .recruitment(JsonFunctions.string(clinicalTrialsEntryObject, "recruitment"))
                            .therapy(retrieveClinicalTrialsTherapies(clinicalTrialsEntryObject.getAsJsonArray("therapies")))
                            .ageGroup(JsonFunctions.stringList(clinicalTrialsEntryObject, "ageGroups"))
                            .gender(JsonFunctions.optionalNullableString(clinicalTrialsEntryObject, "gender"))
                            .variantRequirement(JsonFunctions.string(clinicalTrialsEntryObject, "variantRequirements"))
                            .sponsors(JsonFunctions.optionalNullableString(clinicalTrialsEntryObject, "sponsors"))
                            .updateDate(JsonFunctions.string(clinicalTrialsEntryObject, "updateDate"))
                            .indication(retrieveClinicalTrialsIndications(clinicalTrialsEntryObject.getAsJsonArray("indications")))
                            .variantRequirementDetail(retrieveClinicalTrialsVariantRequirementDetails(clinicalTrialsEntryObject.getAsJsonArray(
                                    "variantRequirementDetails")))
                            .clinicalTrialLocation(retrieveClinicalTrialsLocations(clinicalTrialsEntryObject.getAsJsonArray(
                                    "clinicalTrialLocations")))
                            .build());
                }
                reader.close();
            }
        }
        LOGGER.info("Finished reading clinical trials");

        return clinicalTrials;
    }

    @NotNull
    private static List<ClinicalTrialVariantRequirementDetail> retrieveClinicalTrialsVariantRequirementDetails(
            @NotNull JsonArray jsonArray) {
        List<ClinicalTrialVariantRequirementDetail> variantRequirementDetails = Lists.newArrayList();
        JsonDatamodelChecker clinicalTrailVariantRequirementDetailChecker =
                ClinicalTrialDataModelChecker.clinicalTrialVariantRequirementDetailsObjectChecker();

        for (JsonElement variantRequirementDetail : jsonArray) {
            JsonObject variantRequirementDetailObject = variantRequirementDetail.getAsJsonObject();
            clinicalTrailVariantRequirementDetailChecker.check(variantRequirementDetailObject);

            variantRequirementDetails.add(ImmutableClinicalTrialVariantRequirementDetail.builder()
                    .molecularProfile(retrieveClinicalTrialsMolecularProfile(variantRequirementDetailObject.getAsJsonObject("molecularProfile")))
                    .requirementType(JsonFunctions.string(variantRequirementDetailObject, "requirementType"))
                    .build());
        }
        return variantRequirementDetails;
    }

    @NotNull
    private static MolecularProfileInfo retrieveClinicalTrialsMolecularProfile(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker clinicalTrailMolecularProfileChecker =
                ClinicalTrialDataModelChecker.clinicalTrialMolecularProfileObjectChecker();
        clinicalTrailMolecularProfileChecker.check(jsonObject);

        return ImmutableMolecularProfileInfo.builder()
                .id(JsonFunctions.nullableString(jsonObject, "id"))
                .profileName(JsonFunctions.string(jsonObject, "profileName"))
                .build();
    }

    @NotNull
    private static List<TherapyInfo> retrieveClinicalTrialsTherapies(@NotNull JsonArray jsonArray) {
        List<TherapyInfo> therapies = Lists.newArrayList();
        JsonDatamodelChecker clinicalTrailTherapiesChecker = ClinicalTrialDataModelChecker.clinicalTrialTherapiesObjectChecker();

        for (JsonElement therapy : jsonArray) {
            JsonObject therapyObject = therapy.getAsJsonObject();
            clinicalTrailTherapiesChecker.check(therapyObject);

            therapies.add(ImmutableTherapyInfo.builder()
                    .id(JsonFunctions.string(therapyObject, "id"))
                    .therapyName(JsonFunctions.string(therapyObject, "therapyName"))
                    .synonyms(JsonFunctions.optionalNullableString(therapyObject, "synonyms"))
                    .build());
        }
        return therapies;
    }

    @NotNull
    private static List<IndicationInfo> retrieveClinicalTrialsIndications(@NotNull JsonArray jsonArray) {
        List<IndicationInfo> clinicalTrialsIndications = Lists.newArrayList();
        JsonDatamodelChecker clinicalTrailIndicationsChecker = ClinicalTrialDataModelChecker.clinicalTrialIndicationsObjectChecker();

        for (JsonElement clinicalTrialIndication : jsonArray) {
            JsonObject clinicalTrialIndicationObject = clinicalTrialIndication.getAsJsonObject();
            clinicalTrailIndicationsChecker.check(clinicalTrialIndicationObject);

            clinicalTrialsIndications.add(ImmutableIndicationInfo.builder()
                    .id(JsonFunctions.string(clinicalTrialIndicationObject, "id"))
                    .name(JsonFunctions.string(clinicalTrialIndicationObject, "name"))
                    .source(JsonFunctions.string(clinicalTrialIndicationObject, "source"))
                    .build());
        }
        return clinicalTrialsIndications;
    }

    @NotNull
    private static List<ClinicalTrialLocation> retrieveClinicalTrialsLocations(@NotNull JsonArray jsonArray) {
        List<ClinicalTrialLocation> clinicalTrialLocations = Lists.newArrayList();
        JsonDatamodelChecker clinicalTrailLocationsChecker = ClinicalTrialDataModelChecker.clinicalTrialLocationsObjectChecker();

        for (JsonElement clinicalTrialLocation : jsonArray) {
            JsonObject clinicalTrailLocationObject = clinicalTrialLocation.getAsJsonObject();
            clinicalTrailLocationsChecker.check(clinicalTrailLocationObject);

            clinicalTrialLocations.add(ImmutableClinicalTrialLocation.builder()
                    .nctId(JsonFunctions.string(clinicalTrailLocationObject, "nctId"))
                    .facility(JsonFunctions.optionalNullableString(clinicalTrailLocationObject, "facility"))
                    .city(JsonFunctions.string(clinicalTrailLocationObject, "city"))
                    .country(JsonFunctions.string(clinicalTrailLocationObject, "country"))
                    .status(JsonFunctions.optionalNullableString(clinicalTrailLocationObject, "status"))
                    .state(JsonFunctions.optionalNullableString(clinicalTrailLocationObject, "state"))
                    .zip(JsonFunctions.optionalNullableString(clinicalTrailLocationObject, "zip"))
                    .clinicalTrialContact(retrieveClinicalTrialsContact(clinicalTrailLocationObject.getAsJsonArray("clinicalTrialContacts")))
                    .build());
        }
        return clinicalTrialLocations;
    }

    @NotNull
    private static List<ClinicalTrialContact> retrieveClinicalTrialsContact(@NotNull JsonArray jsonArray) {
        List<ClinicalTrialContact> clinicalTrialContacts = Lists.newArrayList();
        JsonDatamodelChecker clinicalTrailContactChecker = ClinicalTrialDataModelChecker.clinicalTrialContactObjectChecker();

        for (JsonElement clinicalTrialContact : jsonArray) {
            JsonObject clinicalTrialContactObject = clinicalTrialContact.getAsJsonObject();

            clinicalTrailContactChecker.check(clinicalTrialContactObject);
            clinicalTrialContacts.add(ImmutableClinicalTrialContact.builder()
                    .name(JsonFunctions.optionalNullableString(clinicalTrialContactObject, "name"))
                    .email(JsonFunctions.optionalNullableString(clinicalTrialContactObject, "email"))
                    .phone(JsonFunctions.optionalNullableString(clinicalTrialContactObject, "phone"))
                    .phoneExt(JsonFunctions.optionalNullableString(clinicalTrialContactObject, "phoneExt"))
                    .role(JsonFunctions.string(clinicalTrialContactObject, "role"))
                    .build());

        }
        return clinicalTrialContacts;
    }
}
