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
                            .build());
                }
            }
        }
        LOGGER.info("Finished reading molecular profiles");

        return molecularProfiles;
    }

    @NotNull
    public static List<MolecularProfileGeneVariant> extractGeneVariant(@NotNull JsonArray jsonArray) {
        List<MolecularProfileGeneVariant> geneVariants = Lists.newArrayList();

        for (JsonElement geneVariant : jsonArray) {
            JsonObject geneVariantJsonObject = geneVariant.getAsJsonObject();

            geneVariants.add(ImmutableMolecularProfileGeneVariant.builder()
                    .id(JsonFunctions.string(geneVariantJsonObject, "id"))
                    .fullName(JsonFunctions.string(geneVariantJsonObject, "fullName"))
                    .impact(JsonFunctions.nullableString(geneVariantJsonObject, "impact"))
                    .proteinEffect(JsonFunctions.nullableString(geneVariantJsonObject, "proteinEffect"))
                    .build());
        }
        return geneVariants;
    }

    @NotNull
    public static List<MolecularProfileProfileTreatmentApproache> extractProfileTreatmentApproach(@NotNull JsonArray jsonArray) {
        List<MolecularProfileProfileTreatmentApproache> profileTreatmentApproaches = Lists.newArrayList();

        for (JsonElement profileTreatmentApproach : jsonArray) {
            JsonObject profileTreatmentApproachJsonObject = profileTreatmentApproach.getAsJsonObject();

            profileTreatmentApproaches.add(ImmutableMolecularProfileProfileTreatmentApproache.builder()
                    .id(JsonFunctions.string(profileTreatmentApproachJsonObject, "id"))
                    .name(JsonFunctions.string(profileTreatmentApproachJsonObject, "id"))
                    .profileName(JsonFunctions.string(profileTreatmentApproachJsonObject, "id"))
                    .build());
        }
        return profileTreatmentApproaches;
    }

    @NotNull
    public static MolecularProfileComplexMolecularProfileEvidence extractComplexMolecularProfileEvidence(@NotNull JsonObject jsonObject) {
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

        for (JsonElement complexMolecularProfileEvidence : jsonArray) {
            JsonObject complexMolecularProfileEvidenceJsonObject = complexMolecularProfileEvidence.getAsJsonObject();
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
    public static MolecularProfileMolecularProfile extractMolecularProfile(@NotNull JsonObject jsonObject) {
        return ImmutableMolecularProfileMolecularProfile.builder()
                .id(JsonFunctions.string(jsonObject, "id"))
                .profileName(JsonFunctions.string(jsonObject, "profileName"))
                .build();
    }

    @NotNull
    public static MolecularProfileTherapy extractTherapy(@NotNull JsonObject jsonObject) {
        return ImmutableMolecularProfileTherapy.builder()
                .id(JsonFunctions.string(jsonObject, "id"))
                .therapyName(JsonFunctions.string(jsonObject, "therapyName"))
                .synonyms(JsonFunctions.nullableString(jsonObject, "synonyms"))
                .build();
    }

    @NotNull
    public static MolecularProfileIndication extractIndication(@NotNull JsonObject jsonObject) {
        return ImmutableMolecularProfileIndication.builder()
                .id(JsonFunctions.string(jsonObject, "id"))
                .name(JsonFunctions.string(jsonObject, "name"))
                .source(JsonFunctions.string(jsonObject, "source"))
                .build();
    }

    @NotNull
    public static List<MolecularProfileReferences> extractReference(@NotNull JsonArray jsonArray) {
        List<MolecularProfileReferences> references = Lists.newArrayList();

        for (JsonElement reference : jsonArray) {
            JsonObject referenceJsonObject = reference.getAsJsonObject();

            references.add(ImmutableMolecularProfileReferences.builder()
                    .id(JsonFunctions.string(referenceJsonObject, "id"))
                    .pubMedId(JsonFunctions.string(referenceJsonObject, "id"))
                    .title(JsonFunctions.string(referenceJsonObject, "id"))
                    .url(JsonFunctions.string(referenceJsonObject, "id"))
                    .build());
        }
        return references;
    }

    @NotNull
    public static List<MolecularProfileRelevantTreatmentApproach> extractRelevantTreatmentApproach(@NotNull JsonArray jsonArray) {
        List<MolecularProfileRelevantTreatmentApproach> relevantTreatmentApproaches = Lists.newArrayList();

        for (JsonElement relevantTreatmentApproach : jsonArray) {
            JsonObject relevantTreatmentApproachJsonObject = relevantTreatmentApproach.getAsJsonObject();

            relevantTreatmentApproaches.add(ImmutableMolecularProfileRelevantTreatmentApproach.builder()
                    .id(JsonFunctions.string(relevantTreatmentApproachJsonObject, "id"))
                    .name(JsonFunctions.string(relevantTreatmentApproachJsonObject, "name"))
                    .profileName(JsonFunctions.string(relevantTreatmentApproachJsonObject, "profileName"))
                    .build());
        }
        return relevantTreatmentApproaches;
    }
}
