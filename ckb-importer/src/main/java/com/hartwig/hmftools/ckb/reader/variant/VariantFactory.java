package com.hartwig.hmftools.ckb.reader.variant;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.ParseException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;
import com.hartwig.hmftools.ckb.datamodel.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.datamodel.common.EffectInfo;
import com.hartwig.hmftools.ckb.datamodel.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.datamodel.common.GeneInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableDescriptionInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableEffectInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableEvidenceInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableGeneInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableIndicationInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableMolecularProfileInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableReferenceInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableTherapyInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableTreatmentApproachInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableVariantInfo;
import com.hartwig.hmftools.ckb.datamodel.common.IndicationInfo;
import com.hartwig.hmftools.ckb.datamodel.common.MolecularProfileInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ReferenceInfo;
import com.hartwig.hmftools.ckb.datamodel.common.TherapyInfo;
import com.hartwig.hmftools.ckb.datamodel.common.TreatmentApproachInfo;
import com.hartwig.hmftools.ckb.datamodel.common.VariantInfo;
import com.hartwig.hmftools.ckb.datamodel.variant.ImmutableVariant;
import com.hartwig.hmftools.ckb.datamodel.variant.ImmutableVariantCategoryVariantPath;
import com.hartwig.hmftools.ckb.datamodel.variant.ImmutableVariantPartnerGene;
import com.hartwig.hmftools.ckb.datamodel.variant.ImmutableVariantTranscriptCoordinate;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;
import com.hartwig.hmftools.ckb.datamodel.variant.VariantCategoryVariantPath;
import com.hartwig.hmftools.ckb.datamodel.variant.VariantPartnerGene;
import com.hartwig.hmftools.ckb.datamodel.variant.VariantTranscriptCoordinate;
import com.hartwig.hmftools.ckb.util.DateConverter;
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
    public static List<Variant> readingVariant(@NotNull String variantDir) throws IOException, ParseException {
        LOGGER.info("Start reading variant dir");

        List<Variant> variants = Lists.newArrayList();
        File[] filesVariant = new File(variantDir).listFiles();

        if (filesVariant != null) {
            LOGGER.info("The total files in the variant dir is {}", filesVariant.length);

            for (File variant : filesVariant) {
                JsonParser parser = new JsonParser();
                JsonReader reader = new JsonReader(new FileReader(variant));
                reader.setLenient(true);

                while (reader.peek() != JsonToken.END_DOCUMENT) {
                    JsonObject variantEntryObject = parser.parse(reader).getAsJsonObject();
                    JsonDatamodelChecker variantObjectChecker = VariantDataModelChecker.variantObjectChecker();
                    variantObjectChecker.check(variantEntryObject);

                    variants.add(ImmutableVariant.builder()
                            .id(JsonFunctions.integer(variantEntryObject, "id"))
                            .fullName(JsonFunctions.string(variantEntryObject, "fullName"))
                            .impact(JsonFunctions.nullableString(variantEntryObject, "impact"))
                            .proteinEffect(JsonFunctions.nullableString(variantEntryObject, "proteinEffect"))
                            .description(extractGeneDescription(variantEntryObject.getAsJsonArray("geneVariantDescriptions")))
                            .type(JsonFunctions.nullableString(variantEntryObject, "type"))
                            .gene(extractGene(variantEntryObject.getAsJsonObject("gene")))
                            .variant(JsonFunctions.string(variantEntryObject, "variant"))
                            .createDate(DateConverter.convertDate(JsonFunctions.string(variantEntryObject, "createDate")))
                            .updateDate(DateConverter.convertDate(JsonFunctions.string(variantEntryObject, "updateDate")))
                            .referenceTranscriptCoordinate(
                                    variantEntryObject.has("referenceTranscriptCoordinates") && !variantEntryObject.get(
                                            "referenceTranscriptCoordinates").isJsonNull() ? extractReferenceTranscriptCoordinate(
                                            variantEntryObject.getAsJsonObject("referenceTranscriptCoordinates")) : null)
                            .partnerGene(extractPartnerGene(variantEntryObject.getAsJsonArray("partnerGenes")))
                            .categoryVariantPath(extractCategoryVariantPath(variantEntryObject.getAsJsonArray("categoryVariantPaths")))
                            .evidence(extractEvidence(variantEntryObject.getAsJsonArray("evidence")))
                            .extendedEvidence(extractExtendedEvidence(variantEntryObject.getAsJsonArray("extendedEvidence")))
                            .molecularProfile(extractMolecularProfilesList(variantEntryObject.getAsJsonArray("molecularProfiles")))
                            .allTranscriptCoordinate(extractAllTranscriptCoordinates(variantEntryObject.getAsJsonArray(
                                    "allTranscriptCoordinates")))
                            .memberVariant(extractMemberVariants(variantEntryObject.getAsJsonArray("memberVariants")))
                            .build());
                }
                reader.close();
            }
        }
        LOGGER.info("Finished reading variant dir");

        return variants;
    }

    @NotNull
    private static List<DescriptionInfo> extractGeneDescription(@NotNull JsonArray jsonArray) {
        List<DescriptionInfo> geneVariantDescriptions = Lists.newArrayList();
        JsonDatamodelChecker variantGeneVariantDescriptionObjectChecker = VariantDataModelChecker.geneVariantDescriptionObjectChecker();

        for (JsonElement geneVariantDescription : jsonArray) {
            JsonObject geneVariantDescriptionJsonObject = geneVariantDescription.getAsJsonObject();
            variantGeneVariantDescriptionObjectChecker.check(geneVariantDescriptionJsonObject);

            geneVariantDescriptions.add(ImmutableDescriptionInfo.builder()
                    .description(JsonFunctions.string(geneVariantDescriptionJsonObject, "description"))
                    .reference(extractReferences(geneVariantDescriptionJsonObject.getAsJsonArray("references")))
                    .build());
        }
        return geneVariantDescriptions;
    }

    @NotNull
    private static List<ReferenceInfo> extractReferences(@NotNull JsonArray jsonArray) {
        List<ReferenceInfo> references = Lists.newArrayList();
        JsonDatamodelChecker referenceObjectChecker = VariantDataModelChecker.referenceObjectChecker();

        for (JsonElement reference : jsonArray) {
            JsonObject referenceJsonObject = reference.getAsJsonObject();
            referenceObjectChecker.check(referenceJsonObject);

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
    private static GeneInfo extractGene(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker geneObjectChecker = VariantDataModelChecker.geneObjectChecker();
        geneObjectChecker.check(jsonObject);

        return ImmutableGeneInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .geneSymbol(JsonFunctions.string(jsonObject, "geneSymbol"))
                .term(JsonFunctions.stringList(jsonObject, "terms"))
                .build();
    }

    @NotNull
    private static VariantTranscriptCoordinate extractReferenceTranscriptCoordinate(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker geneObjectChecker = VariantDataModelChecker.referenceTranscriptCoordinateObjectChecker();
        geneObjectChecker.check(jsonObject);

        return ImmutableVariantTranscriptCoordinate.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .transcript(JsonFunctions.string(jsonObject, "transcript"))
                .gDNA(JsonFunctions.string(jsonObject, "gDna"))
                .cDNA(JsonFunctions.string(jsonObject, "cDna"))
                .protein(JsonFunctions.string(jsonObject, "protein"))
                .sourceDB(JsonFunctions.string(jsonObject, "sourceDb"))
                .refGenomeBuild(JsonFunctions.string(jsonObject, "refGenomeBuild"))
                .build();
    }

    @NotNull
    private static List<VariantPartnerGene> extractPartnerGene(@NotNull JsonArray jsonArray) {
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
    private static List<VariantCategoryVariantPath> extractCategoryVariantPath(@NotNull JsonArray jsonArray) {
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
    private static List<VariantInfo> extractVariant(@NotNull JsonArray jsonArray) {
        List<VariantInfo> variants = Lists.newArrayList();
        JsonDatamodelChecker variantObjectChecker = VariantDataModelChecker.variantVariantObjectChecker();

        for (JsonElement variant : jsonArray) {
            JsonObject variantJsonObject = variant.getAsJsonObject();
            variantObjectChecker.check(variantJsonObject);

            variants.add(ImmutableVariantInfo.builder()
                    .id(JsonFunctions.integer(variantJsonObject, "id"))
                    .fullName(JsonFunctions.string(variantJsonObject, "fullName"))
                    .impact(JsonFunctions.string(variantJsonObject, "impact"))
                    .proteinEffect(JsonFunctions.string(variantJsonObject, "proteinEffect"))
                    .build());
        }
        return variants;
    }

    @NotNull
    private static List<EvidenceInfo> extractEvidence(@NotNull JsonArray jsonArray) {
        List<EvidenceInfo> evidences = Lists.newArrayList();
        JsonDatamodelChecker evidenceChecker = VariantDataModelChecker.evidenceObjectChecker();

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
                    .reference(extractReferences(evidenceJsonObject.getAsJsonArray("references")))
                    .ampCapAscoEvidenceLevel(JsonFunctions.string(evidenceJsonObject, "ampCapAscoEvidenceLevel"))
                    .ampCapAscoInferredTier(JsonFunctions.string(evidenceJsonObject, "ampCapAscoInferredTier"))
                    .build());
        }
        return evidences;
    }

    @NotNull
    private static MolecularProfileInfo extractMolecularProfile(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker molecularProfileChecker = VariantDataModelChecker.molecularProfileObjectChecker();
        molecularProfileChecker.check(jsonObject);

        return ImmutableMolecularProfileInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .profileName(JsonFunctions.string(jsonObject, "profileName"))
                .build();
    }

    @NotNull
    private static TherapyInfo extractTherapy(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker therapyChecker = VariantDataModelChecker.therapyChecker();
        therapyChecker.check(jsonObject);

        return ImmutableTherapyInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .therapyName(JsonFunctions.string(jsonObject, "therapyName"))
                .synonyms(JsonFunctions.nullableString(jsonObject, "synonyms"))
                .build();
    }

    @NotNull
    private static IndicationInfo extractIndication(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker indicationChecker = VariantDataModelChecker.indicationChecker();
        indicationChecker.check(jsonObject);

        return ImmutableIndicationInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .name(JsonFunctions.string(jsonObject, "name"))
                .source(JsonFunctions.string(jsonObject, "source"))
                .build();
    }

    @NotNull
    private static List<EvidenceInfo> extractExtendedEvidence(@NotNull JsonArray jsonArray) {
        List<EvidenceInfo> extendedEvidences = Lists.newArrayList();
        JsonDatamodelChecker extendedEvidenceChecker = VariantDataModelChecker.extendedEvidenceObjectChecker();

        for (JsonElement extendedEvidence : jsonArray) {
            JsonObject extendedEvidenceJsonObject = extendedEvidence.getAsJsonObject();
            extendedEvidenceChecker.check(extendedEvidenceJsonObject);

            extendedEvidences.add(ImmutableEvidenceInfo.builder()
                    .id(JsonFunctions.integer(extendedEvidenceJsonObject, "id"))
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
    private static List<MolecularProfileInfo> extractMolecularProfilesList(@NotNull JsonArray jsonArray) {
        List<MolecularProfileInfo> molecularProfiles = Lists.newArrayList();
        JsonDatamodelChecker molecularProfileChecker = VariantDataModelChecker.molecularProfileObjectChecker();

        for (JsonElement molecularProfile : jsonArray) {
            JsonObject molecularProfileJsonObject = molecularProfile.getAsJsonObject();
            molecularProfileChecker.check(molecularProfileJsonObject);
            molecularProfiles.add(ImmutableMolecularProfileInfo.builder()
                    .id(JsonFunctions.integer(molecularProfileJsonObject, "id"))
                    .profileName(JsonFunctions.string(molecularProfileJsonObject, "profileName"))
                    .treatmentApproach(extractProfileTreatmentApproches(molecularProfileJsonObject.getAsJsonArray(
                            "profileTreatmentApproaches")))
                    .build());
        }

        return molecularProfiles;
    }

    @NotNull
    private static List<TreatmentApproachInfo> extractProfileTreatmentApproches(@NotNull JsonArray jsonArray) {
        List<TreatmentApproachInfo> profileTreatmentApproaches = Lists.newArrayList();
        JsonDatamodelChecker profileTreatmentApprochChecker = VariantDataModelChecker.profileTreatmentApprochObjectChecker();

        for (JsonElement profileTreatmentApproch : jsonArray) {
            JsonObject profileTreatmentApprochJsonObject = profileTreatmentApproch.getAsJsonObject();
            profileTreatmentApprochChecker.check(profileTreatmentApprochJsonObject);

            profileTreatmentApproaches.add(ImmutableTreatmentApproachInfo.builder()
                    .id(JsonFunctions.integer(profileTreatmentApprochJsonObject, "id"))
                    .name(JsonFunctions.string(profileTreatmentApprochJsonObject, "name"))
                    .profileName(JsonFunctions.string(profileTreatmentApprochJsonObject, "profileName"))
                    .build());
        }

        return profileTreatmentApproaches;
    }

    @NotNull
    private static List<VariantTranscriptCoordinate> extractAllTranscriptCoordinates(@NotNull JsonArray jsonArray) {
        List<VariantTranscriptCoordinate> allTranscriptCoordinates = Lists.newArrayList();
        JsonDatamodelChecker allTranscriptCoordinateChecker = VariantDataModelChecker.allTranscriptCoordinateObjectChecker();

        for (JsonElement allTranscriptCoordinate : jsonArray) {
            JsonObject allTranscriptCoordinatesJsonObject = allTranscriptCoordinate.getAsJsonObject();
            allTranscriptCoordinateChecker.check(allTranscriptCoordinatesJsonObject);

            allTranscriptCoordinates.add(ImmutableVariantTranscriptCoordinate.builder()
                    .id(JsonFunctions.integer(allTranscriptCoordinatesJsonObject, "id"))
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
    private static List<EffectInfo> extractMemberVariants(@NotNull JsonArray jsonArray) {
        List<EffectInfo> memberVariants = Lists.newArrayList();
        JsonDatamodelChecker memberVariantChecker = VariantDataModelChecker.memberVariantObjectChecker();

        for (JsonElement memberVariant : jsonArray) {
            JsonObject memberVariantObject = memberVariant.getAsJsonObject();
            memberVariantChecker.check(memberVariantObject);

            memberVariants.add(ImmutableEffectInfo.builder()
                    .id(JsonFunctions.integer(memberVariantObject, "id"))
                    .fullName(JsonFunctions.string(memberVariantObject, "fullName"))
                    .impact(JsonFunctions.string(memberVariantObject, "impact"))
                    .proteinEffect(JsonFunctions.string(memberVariantObject, "proteinEffect"))
                    .description(extractGeneDescription(memberVariantObject.getAsJsonArray("geneVariantDescriptions")))
                    .build());
        }
        return memberVariants;
    }
}
