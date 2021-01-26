package com.hartwig.hmftools.ckb.molecularprofile;

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
import com.hartwig.hmftools.ckb.common.ImmutableTreatmentApproach;
import com.hartwig.hmftools.ckb.common.ImmutableVariantInfo;
import com.hartwig.hmftools.ckb.common.IndicationInfo;
import com.hartwig.hmftools.ckb.common.MolecularProfileInfo;
import com.hartwig.hmftools.ckb.common.ReferenceInfo;
import com.hartwig.hmftools.ckb.common.TherapyInfo;
import com.hartwig.hmftools.ckb.common.TreatmentApproach;
import com.hartwig.hmftools.ckb.common.VariantInfo;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;
import com.hartwig.hmftools.common.utils.json.JsonFunctions;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class MolecularprofileFactory {

    private static final Logger LOGGER = LogManager.getLogger(MolecularprofileFactory.class);

    private MolecularprofileFactory() {

    }

    @NotNull
    public static List<MolecularProfile> readingMolecularprofile(@NotNull String molecularprofileDir) throws IOException {
        LOGGER.info("Start reading molecular profiles");

        List<MolecularProfile> molecularProfiles = Lists.newArrayList();
        File[] filesMolecularProfiles = new File(molecularprofileDir).listFiles();

        if (filesMolecularProfiles != null) {
            LOGGER.info("The total files in the molecular profiles dir is {}", filesMolecularProfiles.length);

            for (File molecularProfile : filesMolecularProfiles) {
                JsonParser parser = new JsonParser();
                JsonReader reader = new JsonReader(new FileReader(molecularProfile));
                reader.setLenient(true);

                while (reader.peek() != JsonToken.END_DOCUMENT) {
                    JsonObject molecularProfileEntryObject = parser.parse(reader).getAsJsonObject();
                    JsonDatamodelChecker molecularProfileChecker = MolecularProfileDataModelChecker.molecularProfileObjectChecker();
                    molecularProfileChecker.check(molecularProfileEntryObject);

                    molecularProfiles.add(ImmutableMolecularProfile.builder()
                            .id(JsonFunctions.string(molecularProfileEntryObject, "id"))
                            .profileName(JsonFunctions.string(molecularProfileEntryObject, "profileName"))
                            .geneVariant(extractGeneVariant(molecularProfileEntryObject.getAsJsonArray("geneVariants")))
                            .profileProfileTreatmentApproache(extractProfileTreatmentApproach(molecularProfileEntryObject.getAsJsonArray(
                                    "profileTreatmentApproaches")))
                            .createDate(JsonFunctions.string(molecularProfileEntryObject, "createDate"))
                            .updateDate(JsonFunctions.string(molecularProfileEntryObject, "updateDate"))
                            .complexMolecularProfileEvidence(extractComplexMolecularProfileEvidence(molecularProfileEntryObject.getAsJsonObject(
                                    "complexMolecularProfileEvidence")))
                            .treatmentApproachEvidence(extractTreatmentApproachEvidence(molecularProfileEntryObject.getAsJsonObject(
                                    "treatmentApproachEvidence")))
                            .variantAssociatedClinicalTrial(extractVariantAssociatedClinicalTrials(molecularProfileEntryObject.getAsJsonArray(
                                    "variantAssociatedClinicalTrials")))
                            .variantLevelEvidence(extractVariantlevelEvidence(molecularProfileEntryObject.getAsJsonObject("variantLevelEvidence")))
                            .extendedEvidence(extractExtendedEvidence(molecularProfileEntryObject.getAsJsonObject("extendedEvidence")))
                            .build());
                }
                reader.close();
            }
        }
        LOGGER.info("Finished reading molecular profiles");

        return molecularProfiles;
    }

    @NotNull
    public static List<VariantInfo> extractGeneVariant(@NotNull JsonArray jsonArray) {
        List<VariantInfo> geneVariants = Lists.newArrayList();
        JsonDatamodelChecker geneVariantChecker = MolecularProfileDataModelChecker.molecularProfileGeneVariantObjectChecker();

        for (JsonElement geneVariant : jsonArray) {
            JsonObject geneVariantJsonObject = geneVariant.getAsJsonObject();
            geneVariantChecker.check(geneVariantJsonObject);

            geneVariants.add(ImmutableVariantInfo.builder()
                    .id(JsonFunctions.string(geneVariantJsonObject, "id"))
                    .fullName(JsonFunctions.string(geneVariantJsonObject, "fullName"))
                    .impact(JsonFunctions.nullableString(geneVariantJsonObject, "impact"))
                    .proteinEffect(JsonFunctions.nullableString(geneVariantJsonObject, "proteinEffect"))
                    .build());
        }
        return geneVariants;
    }

    @NotNull
    public static List<TreatmentApproach> extractProfileTreatmentApproach(@NotNull JsonArray jsonArray) {
        List<TreatmentApproach> profileTreatmentApproaches = Lists.newArrayList();
        JsonDatamodelChecker profileTreatmentApproachChecker =
                MolecularProfileDataModelChecker.molecularProfileProfileTreatmentApproacjObjectChecker();

        for (JsonElement profileTreatmentApproach : jsonArray) {
            JsonObject profileTreatmentApproachJsonObject = profileTreatmentApproach.getAsJsonObject();
            profileTreatmentApproachChecker.check(profileTreatmentApproachJsonObject);

            profileTreatmentApproaches.add(ImmutableTreatmentApproach.builder()
                    .id(JsonFunctions.string(profileTreatmentApproachJsonObject, "id"))
                    .name(JsonFunctions.string(profileTreatmentApproachJsonObject, "name"))
                    .profileName(JsonFunctions.string(profileTreatmentApproachJsonObject, "profileName"))
                    .build());
        }
        return profileTreatmentApproaches;
    }

    @NotNull
    public static MolecularProfileComplexMolecularProfileEvidence extractComplexMolecularProfileEvidence(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker complexMolecularProfileEvidenceChecker =
                MolecularProfileDataModelChecker.molecularProfileComplexMolecularProfileEvidence();
        complexMolecularProfileEvidenceChecker.check(jsonObject);

        return ImmutableMolecularProfileComplexMolecularProfileEvidence.builder()
                .totalCount(JsonFunctions.string(jsonObject, "totalCount"))
                .complexMolecularProfileEvidence(extractComplexMolecularProfileEvidenceList(jsonObject.getAsJsonArray(
                        "complexMolecularProfileEvidence")))
                .build();
    }

    @NotNull
    public static List<MolecularProfileComplexMolecularProfileEvidenceList> extractComplexMolecularProfileEvidenceList(
            @NotNull JsonArray jsonArray) {
        List<MolecularProfileComplexMolecularProfileEvidenceList> complexMolecularProfileEvidenceList = Lists.newArrayList();
        JsonDatamodelChecker complexMolecularProfileEvidenceListChecker =
                MolecularProfileDataModelChecker.molecularProfileComplexMolecularProfileEvidenceList();

        for (JsonElement complexMolecularProfileEvidence : jsonArray) {
            JsonObject complexMolecularProfileEvidenceJsonObject = complexMolecularProfileEvidence.getAsJsonObject();
            complexMolecularProfileEvidenceListChecker.check(complexMolecularProfileEvidenceJsonObject);

            complexMolecularProfileEvidenceList.add(ImmutableMolecularProfileComplexMolecularProfileEvidenceList.builder()
                    .id(JsonFunctions.string(complexMolecularProfileEvidenceJsonObject, "id"))
                    .approvalStatus(JsonFunctions.string(complexMolecularProfileEvidenceJsonObject, "approvalStatus"))
                    .evidenceType(JsonFunctions.string(complexMolecularProfileEvidenceJsonObject, "evidenceType"))
                    .efficacyEvidence(JsonFunctions.string(complexMolecularProfileEvidenceJsonObject, "efficacyEvidence"))
                    .molecularProfile(extractMolecularProfile(complexMolecularProfileEvidenceJsonObject.getAsJsonObject("molecularProfile")))
                    .therapy(extractTherapy(complexMolecularProfileEvidenceJsonObject.getAsJsonObject("therapy")))
                    .indication(extractIndication(complexMolecularProfileEvidenceJsonObject.getAsJsonObject("indication")))
                    .responseType(JsonFunctions.string(complexMolecularProfileEvidenceJsonObject, "responseType"))
                    .reference(extractReference(complexMolecularProfileEvidenceJsonObject.getAsJsonArray("references")))
                    .ampCapAscoEvidenceLevel(JsonFunctions.string(complexMolecularProfileEvidenceJsonObject, "ampCapAscoEvidenceLevel"))
                    .ampCapAscoInferredTier(JsonFunctions.string(complexMolecularProfileEvidenceJsonObject, "ampCapAscoInferredTier"))
                    .relevantTreatmentApproach(extractRelevantTreatmentApproach(complexMolecularProfileEvidenceJsonObject.getAsJsonArray(
                            "relevantTreatmentApproaches")))
                    .build());
        }
        return complexMolecularProfileEvidenceList;
    }

    @NotNull
    public static MolecularProfileInfo extractMolecularProfile(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker molecularProfileChecker = MolecularProfileDataModelChecker.molecularProfileMolecularprofile();
        molecularProfileChecker.check(jsonObject);

        return ImmutableMolecularProfileInfo.builder()
                .id(JsonFunctions.string(jsonObject, "id"))
                .profileName(JsonFunctions.string(jsonObject, "profileName"))
                .build();
    }

    @NotNull
    public static TherapyInfo extractTherapy(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker therapyChecker = MolecularProfileDataModelChecker.molecularProfileTherapy();
        therapyChecker.check(jsonObject);

        return ImmutableTherapyInfo.builder()
                .id(JsonFunctions.string(jsonObject, "id"))
                .therapyName(JsonFunctions.string(jsonObject, "therapyName"))
                .synonyms(JsonFunctions.nullableString(jsonObject, "synonyms"))
                .build();
    }

    @NotNull
    public static IndicationInfo extractIndication(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker indicationChecker = MolecularProfileDataModelChecker.molecularProfileIndication();
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
        JsonDatamodelChecker referenceChecker = MolecularProfileDataModelChecker.molecularProfileReference();

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
    public static List<TreatmentApproach> extractRelevantTreatmentApproach(@NotNull JsonArray jsonArray) {
        List<TreatmentApproach> relevantTreatmentApproaches = Lists.newArrayList();
        JsonDatamodelChecker relevantTreatmentApproachChecker =
                MolecularProfileDataModelChecker.molecularProfileRelevantTreatmentApproach();

        for (JsonElement relevantTreatmentApproach : jsonArray) {
            JsonObject relevantTreatmentApproachJsonObject = relevantTreatmentApproach.getAsJsonObject();
            relevantTreatmentApproachChecker.check(relevantTreatmentApproachJsonObject);

            relevantTreatmentApproaches.add(ImmutableTreatmentApproach.builder()
                    .id(JsonFunctions.string(relevantTreatmentApproachJsonObject, "id"))
                    .name(JsonFunctions.string(relevantTreatmentApproachJsonObject, "name"))
                    .profileName(JsonFunctions.string(relevantTreatmentApproachJsonObject, "profileName"))
                    .build());
        }
        return relevantTreatmentApproaches;
    }

    @NotNull
    public static MolecularProfileTreatmentApproachEvidence extractTreatmentApproachEvidence(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker treatmentApproachEvidenceeChecker =
                MolecularProfileDataModelChecker.molecularProfileTreatmentApproachEvidence();
        treatmentApproachEvidenceeChecker.check(jsonObject);

        return ImmutableMolecularProfileTreatmentApproachEvidence.builder()
                .totalCount(JsonFunctions.string(jsonObject, "totalCount"))
                .treatmentApproachEvidenceList(extractTreatmentApproachEvidenceList(jsonObject.getAsJsonArray("treatmentApproachEvidence")))
                .build();
    }

    @NotNull
    public static List<MolecularProfileTreatmentApproachEvidenceList> extractTreatmentApproachEvidenceList(@NotNull JsonArray jsonArray) {
        List<MolecularProfileTreatmentApproachEvidenceList> treatmentApproachEvidenceList = Lists.newArrayList();
        JsonDatamodelChecker treatmentApproachEvidenceListChecker =
                MolecularProfileDataModelChecker.molecularProfileTreatmentApproachEvidenceList();

        for (JsonElement treatmentApproachEvidence : jsonArray) {
            JsonObject treatmentApproachEvidenceJsonObject = treatmentApproachEvidence.getAsJsonObject();
            treatmentApproachEvidenceListChecker.check(treatmentApproachEvidenceJsonObject);

            treatmentApproachEvidenceList.add(ImmutableMolecularProfileTreatmentApproachEvidenceList.builder()
                    .id(JsonFunctions.string(treatmentApproachEvidenceJsonObject, "id"))
                    .approvalStatus(JsonFunctions.string(treatmentApproachEvidenceJsonObject, "approvalStatus"))
                    .evidenceType(JsonFunctions.string(treatmentApproachEvidenceJsonObject, "evidenceType"))
                    .efficacyEvidence(JsonFunctions.string(treatmentApproachEvidenceJsonObject, "efficacyEvidence"))
                    .molecularProfile(extractMolecularProfile(treatmentApproachEvidenceJsonObject.getAsJsonObject("molecularProfile")))
                    .therapy(extractTherapy(treatmentApproachEvidenceJsonObject.getAsJsonObject("therapy")))
                    .indication(extractIndication(treatmentApproachEvidenceJsonObject.getAsJsonObject("indication")))
                    .responseType(JsonFunctions.string(treatmentApproachEvidenceJsonObject, "responseType"))
                    .reference(extractReference(treatmentApproachEvidenceJsonObject.getAsJsonArray("references")))
                    .ampCapAscoEvidenceLevel(JsonFunctions.string(treatmentApproachEvidenceJsonObject, "ampCapAscoEvidenceLevel"))
                    .ampCapAscoInferredTier(JsonFunctions.string(treatmentApproachEvidenceJsonObject, "ampCapAscoInferredTier"))
                    .relevantTreatmentApproach(extractRelevantTreatmentApproach(treatmentApproachEvidenceJsonObject.getAsJsonArray(
                            "relevantTreatmentApproaches")))
                    .build());
        }
        return treatmentApproachEvidenceList;
    }

    @NotNull
    public static List<MolecularProfileVariantAssociatedClinicalTrials> extractVariantAssociatedClinicalTrials(
            @NotNull JsonArray jsonArray) {
        List<MolecularProfileVariantAssociatedClinicalTrials> variantAssociatedClinicalTrials = Lists.newArrayList();
        JsonDatamodelChecker variantAssociatedClinicalTrialChecker = MolecularProfileDataModelChecker.variantAssociatedClinicalTrial();

        for (JsonElement variantAssociatedClinicalTrial : jsonArray) {
            JsonObject variantAssociatedClinicalTrialJsonObject = variantAssociatedClinicalTrial.getAsJsonObject();
            variantAssociatedClinicalTrialChecker.check(variantAssociatedClinicalTrialJsonObject);

            variantAssociatedClinicalTrials.add(ImmutableMolecularProfileVariantAssociatedClinicalTrials.builder()
                    .nctId(JsonFunctions.string(variantAssociatedClinicalTrialJsonObject, "nctId"))
                    .title(JsonFunctions.string(variantAssociatedClinicalTrialJsonObject, "title"))
                    .phase(JsonFunctions.string(variantAssociatedClinicalTrialJsonObject, "phase"))
                    .recruitment(JsonFunctions.string(variantAssociatedClinicalTrialJsonObject, "recruitment"))
                    .therapy(extractTherapyList(variantAssociatedClinicalTrialJsonObject.getAsJsonArray("therapies")))
                    .build());
        }
        return variantAssociatedClinicalTrials;
    }

    @NotNull
    public static List<TherapyInfo> extractTherapyList(@NotNull JsonArray jsonArray) {
        List<TherapyInfo> therapies = Lists.newArrayList();
        JsonDatamodelChecker therapyChecker = MolecularProfileDataModelChecker.molecularProfileTherapy();

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
    public static MolecularProfileVariantLevelEvidence extractVariantlevelEvidence(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker variantLevelEvidenceChecker =
                MolecularProfileDataModelChecker.molecularProfilevariantLevelEvidence();
        variantLevelEvidenceChecker.check(jsonObject);

        return ImmutableMolecularProfileVariantLevelEvidence.builder()
                .totalCount(JsonFunctions.string(jsonObject, "totalCount"))
                .variantLevelEvidenceList(extractVariantlevelEvidenceList(jsonObject.getAsJsonArray("variantLevelEvidences")))
                .build();
    }

    @NotNull
    public static List<MolecularProfileVariantLevelEvidenceList> extractVariantlevelEvidenceList(@NotNull JsonArray jsonArray) {
        List<MolecularProfileVariantLevelEvidenceList> variantlevelEvidenceList = Lists.newArrayList();
        JsonDatamodelChecker variantlevelEvidenceListChecker =
                MolecularProfileDataModelChecker.molecularProfilevariantLevelEvidenceList();

        for (JsonElement variantLevelEvidence : jsonArray) {
            JsonObject variantLevelEvidenceJsonObject = variantLevelEvidence.getAsJsonObject();
            variantlevelEvidenceListChecker.check(variantLevelEvidenceJsonObject);

            variantlevelEvidenceList.add(ImmutableMolecularProfileVariantLevelEvidenceList.builder()
                    .id(JsonFunctions.string(variantLevelEvidenceJsonObject, "id"))
                    .approvalStatus(JsonFunctions.string(variantLevelEvidenceJsonObject, "approvalStatus"))
                    .evidenceType(JsonFunctions.string(variantLevelEvidenceJsonObject, "evidenceType"))
                    .efficacyEvidence(JsonFunctions.string(variantLevelEvidenceJsonObject, "efficacyEvidence"))
                    .molecularProfile(extractMolecularProfile(variantLevelEvidenceJsonObject.getAsJsonObject("molecularProfile")))
                    .therapy(extractTherapy(variantLevelEvidenceJsonObject.getAsJsonObject("therapy")))
                    .indication(extractIndication(variantLevelEvidenceJsonObject.getAsJsonObject("indication")))
                    .responseType(JsonFunctions.string(variantLevelEvidenceJsonObject, "responseType"))
                    .reference(extractReference(variantLevelEvidenceJsonObject.getAsJsonArray("references")))
                    .ampCapAscoEvidenceLevel(JsonFunctions.string(variantLevelEvidenceJsonObject, "ampCapAscoEvidenceLevel"))
                    .ampCapAscoInferredTier(JsonFunctions.string(variantLevelEvidenceJsonObject, "ampCapAscoInferredTier"))
                    .relevantTreatmentApproach(extractRelevantTreatmentApproach(variantLevelEvidenceJsonObject.getAsJsonArray(
                            "relevantTreatmentApproaches")))
                    .build());
        }
        return variantlevelEvidenceList;
    }

    @NotNull
    public static MolecularProfileExtendedEvidence extractExtendedEvidence(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker extendedEvidenceChecker =
                MolecularProfileDataModelChecker.extendedEvidenceEvidence();
        extendedEvidenceChecker.check(jsonObject);

        return ImmutableMolecularProfileExtendedEvidence.builder()
                .totalCount(JsonFunctions.string(jsonObject, "totalCount"))
                .extendedEvidenceList(extractExtendedEvidenceList(jsonObject.getAsJsonArray("extendedEvidence")))
                .build();
    }

    @NotNull
    public static List<MolecularProfileExtendedEvidenceList> extractExtendedEvidenceList(@NotNull JsonArray jsonArray) {
        List<MolecularProfileExtendedEvidenceList> extendedEvidenceList = Lists.newArrayList();
        JsonDatamodelChecker extendedEvidenceListChecker =
                MolecularProfileDataModelChecker.extendedEvidenceList();

        for (JsonElement extendedEvidence : jsonArray) {
            JsonObject extendedEvidenceJsonObject = extendedEvidence.getAsJsonObject();
            extendedEvidenceListChecker.check(extendedEvidenceJsonObject);

            extendedEvidenceList.add(ImmutableMolecularProfileExtendedEvidenceList.builder()
                    .id(JsonFunctions.string(extendedEvidenceJsonObject, "id"))
                    .approvalStatus(JsonFunctions.string(extendedEvidenceJsonObject, "approvalStatus"))
                    .evidenceType(JsonFunctions.string(extendedEvidenceJsonObject, "evidenceType"))
                    .efficacyEvidence(JsonFunctions.string(extendedEvidenceJsonObject, "efficacyEvidence"))
                    .molecularProfile(extractMolecularProfile(extendedEvidenceJsonObject.getAsJsonObject("molecularProfile")))
                    .therapy(extractTherapy(extendedEvidenceJsonObject.getAsJsonObject("therapy")))
                    .indication(extractIndication(extendedEvidenceJsonObject.getAsJsonObject("indication")))
                    .responseType(JsonFunctions.string(extendedEvidenceJsonObject, "responseType"))
                    .reference(extractReference(extendedEvidenceJsonObject.getAsJsonArray("references")))
                    .ampCapAscoEvidenceLevel(JsonFunctions.string(extendedEvidenceJsonObject, "ampCapAscoEvidenceLevel"))
                    .ampCapAscoInferredTier(JsonFunctions.string(extendedEvidenceJsonObject, "ampCapAscoInferredTier"))

                    .build());
        }
        return extendedEvidenceList;
    }
}
