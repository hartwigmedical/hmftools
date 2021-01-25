package com.hartwig.hmftools.ckb.variant;

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
import com.hartwig.hmftools.ckb.common.ImmutableTherapyInfo;
import com.hartwig.hmftools.ckb.common.TherapyInfo;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;
import com.hartwig.hmftools.common.utils.json.JsonFunctions;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class VariantFactory {

    private static final Logger LOGGER = LogManager.getLogger(VariantFactory.class);

    private VariantFactory() {

    }

    @NotNull
    public static List<Variant> readingVariant(@NotNull String variantDir) throws IOException {
        LOGGER.info("Start reading variant");

        List<Variant> variants = Lists.newArrayList();
        File[] filesVariant = new File(variantDir).listFiles();

        if (filesVariant != null) {
            LOGGER.info("The total files in the varinat dir is {}", filesVariant.length);

            for (File variant : filesVariant) {
                JsonParser parser = new JsonParser();
                JsonReader reader = new JsonReader(new FileReader(variant));
                reader.setLenient(true);

                while (reader.peek() != JsonToken.END_DOCUMENT) {
                    JsonObject variantEntryObject = parser.parse(reader).getAsJsonObject();
                    JsonDatamodelChecker variantObjectChecker = VariantDataModelChecker.variantObjectChecker();
                    variantObjectChecker.check(variantEntryObject);

                    variants.add(ImmutableVariant.builder()
                            .id(JsonFunctions.string(variantEntryObject, "id"))
                            .fullName(JsonFunctions.string(variantEntryObject, "fullName"))
                            .impact(JsonFunctions.nullableString(variantEntryObject, "impact"))
                            .proteinEffect(JsonFunctions.nullableString(variantEntryObject, "proteinEffect"))
                            .geneVariantDescription(extractGeneDescription(variantEntryObject.getAsJsonArray("geneVariantDescriptions")))
                            .type(JsonFunctions.nullableString(variantEntryObject, "type"))
                            .gene(extractGene(variantEntryObject.getAsJsonObject("gene")))
                            .variant(JsonFunctions.string(variantEntryObject, "variant"))
                            .createDate(JsonFunctions.string(variantEntryObject, "createDate"))
                            .updateDate(JsonFunctions.string(variantEntryObject, "updateDate"))
                            .referenceTranscriptCoordinates(
                                    variantEntryObject.has("referenceTranscriptCoordinates") && !variantEntryObject.get(
                                            "referenceTranscriptCoordinates").isJsonNull() ? extractReferenceTranscriptCoordinate(
                                            variantEntryObject.getAsJsonObject("referenceTranscriptCoordinates")) : null)
                            .partnerGene(extractPartnerGene(variantEntryObject.getAsJsonArray("partnerGenes")))
                            .categoryVariantPath(extractCategoryVariantPath(variantEntryObject.getAsJsonArray("categoryVariantPaths")))
                            .evidence(extractEvidence(variantEntryObject.getAsJsonArray("evidence")))
                            .extendedEvidence(extractExtendedEvidence(variantEntryObject.getAsJsonArray("extendedEvidence")))
                            .molecularProfile(extarctMolecularProfilesList(variantEntryObject.getAsJsonArray("molecularProfiles")))
                            .allTranscriptCoordinate(extractAllTranscriptCoordinates(variantEntryObject.getAsJsonArray(
                                    "allTranscriptCoordinates")))
                            .memberVariant(extractMemberVariants(variantEntryObject.getAsJsonArray("memberVariants")))
                            .build());
                }
                reader.close();
            }
        }
        LOGGER.info("Finished reading variant");

        return variants;
    }

    @NotNull
    public static List<VariantGeneVariantDescription> extractGeneDescription(@NotNull JsonArray jsonArray) {

        List<VariantGeneVariantDescription> geneVariantDescriptions = Lists.newArrayList();
        JsonDatamodelChecker variantGeneVariantDescriptionObjectChecker = VariantDataModelChecker.geneVariantDescriptionObjectChecker();

        for (JsonElement geneVariantDescription : jsonArray) {
            JsonObject geneVariantDescriptionJsonObject = geneVariantDescription.getAsJsonObject();
            variantGeneVariantDescriptionObjectChecker.check(geneVariantDescriptionJsonObject);

            geneVariantDescriptions.add(ImmutableVariantGeneVariantDescription.builder()
                    .description(JsonFunctions.string(geneVariantDescriptionJsonObject, "description"))
                    .reference(extractReferences(geneVariantDescriptionJsonObject.getAsJsonArray("references")))
                    .build());
        }

        return geneVariantDescriptions;
    }

    @NotNull
    public static List<VarinatReference> extractReferences(@NotNull JsonArray jsonArray) {
        List<VarinatReference> references = Lists.newArrayList();
        JsonDatamodelChecker referenceObjectChecker = VariantDataModelChecker.referenceObjectChecker();

        for (JsonElement reference : jsonArray) {
            JsonObject referenceJsonObject = reference.getAsJsonObject();
            referenceObjectChecker.check(referenceJsonObject);

            references.add(ImmutableVarinatReference.builder()
                    .id(JsonFunctions.string(referenceJsonObject, "id"))
                    .pubMedId(JsonFunctions.nullableString(referenceJsonObject, "pubMedId"))
                    .title(JsonFunctions.nullableString(referenceJsonObject, "title"))
                    .url(JsonFunctions.nullableString(referenceJsonObject, "url"))
                    .build());
        }

        return references;
    }

    @NotNull
    public static VariantGene extractGene(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker geneObjectChecker = VariantDataModelChecker.geneObjectChecker();
        geneObjectChecker.check(jsonObject);

        return ImmutableVariantGene.builder()
                .id(JsonFunctions.string(jsonObject, "id"))
                .geneSymbol(JsonFunctions.string(jsonObject, "geneSymbol"))
                .terms(JsonFunctions.stringList(jsonObject, "terms"))
                .build();
    }

    @NotNull
    public static VariantReferenceTranscriptCoordinate extractReferenceTranscriptCoordinate(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker geneObjectChecker = VariantDataModelChecker.referenceTranscriptCoordinateObjectChecker();
        geneObjectChecker.check(jsonObject);

        return ImmutableVariantReferenceTranscriptCoordinate.builder()
                .id(JsonFunctions.string(jsonObject, "id"))
                .transcript(JsonFunctions.string(jsonObject, "transcript"))
                .gDNA(JsonFunctions.string(jsonObject, "gDna"))
                .cDNA(JsonFunctions.string(jsonObject, "cDna"))
                .protein(JsonFunctions.string(jsonObject, "protein"))
                .sourceDB(JsonFunctions.string(jsonObject, "sourceDb"))
                .refGenomeBuild(JsonFunctions.string(jsonObject, "refGenomeBuild"))
                .build();
    }

    @NotNull
    public static List<VariantPartnerGene> extractPartnerGene(@NotNull JsonArray jsonArray) {
        List<VariantPartnerGene> partnerGenes = Lists.newArrayList();
        JsonDatamodelChecker partnerGeneObjectChecker = VariantDataModelChecker.partnerGeneObjectChecker();

        for (JsonElement partnerGene : jsonArray) {
            JsonObject partnerGenePathJsonObject = partnerGene.getAsJsonObject();
            partnerGeneObjectChecker.check(partnerGenePathJsonObject);

            partnerGenes.add(ImmutableVariantPartnerGene.builder()
                    .gene(extractGene(partnerGenePathJsonObject.getAsJsonObject("gene")))
                    .build());
        }

        return partnerGenes;
    }

    @NotNull
    public static List<VariantCategoryVariantPath> extractCategoryVariantPath(@NotNull JsonArray jsonArray) {
        List<VariantCategoryVariantPath> categoryVariantPaths = Lists.newArrayList();
        JsonDatamodelChecker categoryVariantPathObjectChecker = VariantDataModelChecker.categoryVariantPathObjectChecker();

        for (JsonElement categoryVariantPath : jsonArray) {
            JsonObject categoryVariantPathJsonObject = categoryVariantPath.getAsJsonObject();
            categoryVariantPathObjectChecker.check(categoryVariantPathJsonObject);

            categoryVariantPaths.add(ImmutableVariantCategoryVariantPath.builder()
                    .variantPath(JsonFunctions.string(categoryVariantPathJsonObject, "variantPath"))
                    .variant(extractVariant(categoryVariantPathJsonObject.getAsJsonArray("variants")))
                    .build());
        }

        return categoryVariantPaths;
    }

    @NotNull
    public static List<VariantVariant> extractVariant(@NotNull JsonArray jsonArray) {
        List<VariantVariant> variants = Lists.newArrayList();
        JsonDatamodelChecker variantObjectChecker = VariantDataModelChecker.variantVariantObjectChecker();

        for (JsonElement variant : jsonArray) {
            JsonObject variantJsonObject = variant.getAsJsonObject();
            variantObjectChecker.check(variantJsonObject);

            variants.add(ImmutableVariantVariant.builder()
                    .id(JsonFunctions.string(variantJsonObject, "id"))
                    .fullName(JsonFunctions.string(variantJsonObject, "fullName"))
                    .impact(JsonFunctions.string(variantJsonObject, "impact"))
                    .proteinEffect(JsonFunctions.string(variantJsonObject, "proteinEffect"))
                    .build());
        }

        return variants;
    }

    @NotNull
    public static List<VariantEvidence> extractEvidence(@NotNull JsonArray jsonArray) {
        List<VariantEvidence> evidences = Lists.newArrayList();
        JsonDatamodelChecker evidenceChecker = VariantDataModelChecker.evidenceObjectChecker();

        for (JsonElement evidence : jsonArray) {
            JsonObject evidenceJsonObject = evidence.getAsJsonObject();
            evidenceChecker.check(evidenceJsonObject);

            evidences.add(ImmutableVariantEvidence.builder()
                    .id(JsonFunctions.string(evidenceJsonObject, "id"))
                    .approvalStatus(JsonFunctions.string(evidenceJsonObject, "approvalStatus"))
                    .evidenceType(JsonFunctions.string(evidenceJsonObject, "evidenceType"))
                    .efficacyEvidence(JsonFunctions.string(evidenceJsonObject, "efficacyEvidence"))
                    .molecularProfile(extractMolecularProfile(evidenceJsonObject.getAsJsonObject("molecularProfile")))
                    .therapy(extractTherapy(evidenceJsonObject.getAsJsonObject("therapy")))
                    .indication(extractIndication(evidenceJsonObject.getAsJsonObject("indication")))
                    .resonseType(JsonFunctions.string(evidenceJsonObject, "responseType"))
                    .reference(extractReferences(evidenceJsonObject.getAsJsonArray("references")))
                    .ampCapAscoEvidenceLevel(JsonFunctions.string(evidenceJsonObject, "ampCapAscoEvidenceLevel"))
                    .ampCapAscoInferredTier(JsonFunctions.string(evidenceJsonObject, "ampCapAscoInferredTier"))
                    .build());
        }
        return evidences;
    }

    @NotNull
    public static VariantMolecularProfile extractMolecularProfile(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker molecularProfileChecker = VariantDataModelChecker.molecularProfileObjectChecker();
        molecularProfileChecker.check(jsonObject);

        return ImmutableVariantMolecularProfile.builder()
                .id(JsonFunctions.string(jsonObject, "id"))
                .profileName(JsonFunctions.string(jsonObject, "profileName"))
                .profileTreatmentApproach(null)
                .build();
    }

    @NotNull
    public static TherapyInfo extractTherapy(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker therapyChecker = VariantDataModelChecker.therapyChecker();
        therapyChecker.check(jsonObject);

        return ImmutableTherapyInfo.builder()
                .id(JsonFunctions.string(jsonObject, "id"))
                .therapyName(JsonFunctions.string(jsonObject, "therapyName"))
                .synonyms(JsonFunctions.nullableString(jsonObject, "synonyms"))
                .build();
    }

    @NotNull
    public static VariantIndication extractIndication(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker indicationChecker = VariantDataModelChecker.indicationChecker();
        indicationChecker.check(jsonObject);

        return ImmutableVariantIndication.builder()
                .id(JsonFunctions.string(jsonObject, "id"))
                .name(JsonFunctions.string(jsonObject, "name"))
                .source(JsonFunctions.string(jsonObject, "source"))
                .build();
    }

    @NotNull
    public static List<VariantExtendedEvidence> extractExtendedEvidence(@NotNull JsonArray jsonArray) {
        List<VariantExtendedEvidence> extendedEvidences = Lists.newArrayList();
        JsonDatamodelChecker extendedEvidenceChecker = VariantDataModelChecker.extendedEvidenceObjectChecker();

        for (JsonElement extendedEvidence : jsonArray) {
            JsonObject extendedEvidenceJsonObject = extendedEvidence.getAsJsonObject();
            extendedEvidenceChecker.check(extendedEvidenceJsonObject);

            extendedEvidences.add(ImmutableVariantExtendedEvidence.builder()
                    .id(JsonFunctions.string(extendedEvidenceJsonObject, "id"))
                    .approvalStatus(JsonFunctions.string(extendedEvidenceJsonObject, "approvalStatus"))
                    .evidenceType(JsonFunctions.string(extendedEvidenceJsonObject, "evidenceType"))
                    .efficacyEvidence(JsonFunctions.string(extendedEvidenceJsonObject, "efficacyEvidence"))
                    .molecularProfile(extractMolecularProfile(extendedEvidenceJsonObject.getAsJsonObject("molecularProfile")))
                    .therapy(extractTherapy(extendedEvidenceJsonObject.getAsJsonObject("therapy")))
                    .indication(extractIndication(extendedEvidenceJsonObject.getAsJsonObject("indication")))
                    .responseType(JsonFunctions.string(extendedEvidenceJsonObject, "responseType"))
                    .reference(extractReferences(extendedEvidenceJsonObject.getAsJsonArray("references")))
                    .ampCapAscoEvidenceLevel(JsonFunctions.string(extendedEvidenceJsonObject, "ampCapAscoEvidenceLevel"))
                    .ampCapAscoInferredTier(JsonFunctions.string(extendedEvidenceJsonObject, "ampCapAscoInferredTier"))
                    .build());
        }
        return extendedEvidences;
    }

    @NotNull
    public static List<VariantMolecularProfile> extarctMolecularProfilesList(@NotNull JsonArray jsonArray) {
        List<VariantMolecularProfile> molecularProfiles = Lists.newArrayList();
        JsonDatamodelChecker molecularProfileChecker = VariantDataModelChecker.molecularProfileObjectChecker();

        for (JsonElement molecularProfile : jsonArray) {
            JsonObject molecularProfileJsonObject = molecularProfile.getAsJsonObject();
            molecularProfileChecker.check(molecularProfileJsonObject);
            molecularProfiles.add(ImmutableVariantMolecularProfile.builder()
                    .id(JsonFunctions.string(molecularProfileJsonObject, "id"))
                    .profileName(JsonFunctions.string(molecularProfileJsonObject, "profileName"))
                    .profileTreatmentApproach(extractProfileTreatmentApproches(molecularProfileJsonObject.getAsJsonArray(
                            "profileTreatmentApproaches")))
                    .build());
        }

        return molecularProfiles;
    }

    @NotNull
    public static List<VariantProfileTreatmentApproach> extractProfileTreatmentApproches(@NotNull JsonArray jsonArray) {
        List<VariantProfileTreatmentApproach> profileTreatmentApproaches = Lists.newArrayList();
        JsonDatamodelChecker profileTreatmentApprochChecker = VariantDataModelChecker.profileTreatmentApprochObjectChecker();

        for (JsonElement profileTreatmentApproch : jsonArray) {
            JsonObject profileTreatmentApprochJsonObject = profileTreatmentApproch.getAsJsonObject();
            profileTreatmentApprochChecker.check(profileTreatmentApprochJsonObject);

            profileTreatmentApproaches.add(ImmutableVariantProfileTreatmentApproach.builder()
                    .id(JsonFunctions.string(profileTreatmentApprochJsonObject, "id"))
                    .name(JsonFunctions.string(profileTreatmentApprochJsonObject, "name"))
                    .profileName(JsonFunctions.string(profileTreatmentApprochJsonObject, "profileName"))
                    .build());
        }

        return profileTreatmentApproaches;
    }

    @NotNull
    public static List<VariantAllTranscriptCoordinate> extractAllTranscriptCoordinates(@NotNull JsonArray jsonArray) {
        List<VariantAllTranscriptCoordinate> allTranscriptCoordinates = Lists.newArrayList();
        JsonDatamodelChecker allTranscriptCoordinateChecker = VariantDataModelChecker.allTranscriptCoordinateObjectChecker();

        for (JsonElement allTranscriptCoordinate : jsonArray) {
            JsonObject allTranscriptCoordinatesJsonObject = allTranscriptCoordinate.getAsJsonObject();
            allTranscriptCoordinateChecker.check(allTranscriptCoordinatesJsonObject);

            allTranscriptCoordinates.add(ImmutableVariantAllTranscriptCoordinate.builder()
                    .id(JsonFunctions.string(allTranscriptCoordinatesJsonObject, "id"))
                    .transcript(JsonFunctions.string(allTranscriptCoordinatesJsonObject, "transcript"))
                    .gDNA(JsonFunctions.string(allTranscriptCoordinatesJsonObject, "gDna"))
                    .cDNA(JsonFunctions.string(allTranscriptCoordinatesJsonObject, "cDna"))
                    .protein(JsonFunctions.string(allTranscriptCoordinatesJsonObject, "protein"))
                    .sourceDB(JsonFunctions.string(allTranscriptCoordinatesJsonObject, "sourceDb"))
                    .refGenomeBuild(JsonFunctions.string(allTranscriptCoordinatesJsonObject, "refGenomeBuild"))
                    .build());
        }

        return allTranscriptCoordinates;
    }

    @NotNull
    public static List<VariantMemberVariant> extractMemberVariants(@NotNull JsonArray jsonArray) {
        List<VariantMemberVariant> memberVariants = Lists.newArrayList();
        JsonDatamodelChecker memberVariantChecker = VariantDataModelChecker.memberVariantObjectChecker();

        for (JsonElement memberVariant : jsonArray) {
            JsonObject memberVariantObject = memberVariant.getAsJsonObject();
            memberVariantChecker.check(memberVariantObject);

            memberVariants.add(ImmutableVariantMemberVariant.builder()
                    .id(JsonFunctions.string(memberVariantObject, "id"))
                    .fullName(JsonFunctions.string(memberVariantObject, "fullName"))
                    .impact(JsonFunctions.string(memberVariantObject, "impact"))
                    .proteinEffect(JsonFunctions.string(memberVariantObject, "proteinEffect"))
                    .geneVariantDescription(extractGeneDescription(memberVariantObject.getAsJsonArray("geneVariantDescriptions")))
                    .build());
        }
        return memberVariants;
    }
}
