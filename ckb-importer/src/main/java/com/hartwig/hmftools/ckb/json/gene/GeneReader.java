package com.hartwig.hmftools.ckb.json.gene;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.ckb.json.CkbJsonDirectoryReader;
import com.hartwig.hmftools.ckb.json.common.ClinicalTrialInfo;
import com.hartwig.hmftools.ckb.json.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.json.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableClinicalTrialInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableDescriptionInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableEvidenceInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableIndicationInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableMolecularProfileInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableReferenceInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableTherapyInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableTreatmentApproachInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableVariantInfo;
import com.hartwig.hmftools.ckb.json.common.IndicationInfo;
import com.hartwig.hmftools.ckb.json.common.MolecularProfileInfo;
import com.hartwig.hmftools.ckb.json.common.ReferenceInfo;
import com.hartwig.hmftools.ckb.json.common.TherapyInfo;
import com.hartwig.hmftools.ckb.json.common.TreatmentApproachInfo;
import com.hartwig.hmftools.ckb.json.common.VariantInfo;
import com.hartwig.hmftools.ckb.util.DateConverter;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;
import com.hartwig.hmftools.common.utils.json.JsonFunctions;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class GeneReader extends CkbJsonDirectoryReader<JsonGene> {
    private static final Logger LOGGER = LogManager.getLogger(GeneReader.class);
    public GeneReader(@Nullable final Integer maxFilesToRead) {
        super(maxFilesToRead);
    }

    @NotNull
    @Override
    protected JsonGene read(@NotNull final JsonObject object) {
        JsonDatamodelChecker geneChecker = GeneDataModelChecker.geneObjectChecker();
        geneChecker.check(object);

        return ImmutableJsonGene.builder()
                .id(JsonFunctions.integer(object, "id"))
                .geneSymbol(JsonFunctions.string(object, "geneSymbol"))
                .terms(JsonFunctions.stringList(object, "terms"))
                .entrezId(JsonFunctions.nullableString(object, "entrezId"))
                .synonyms(JsonFunctions.stringList(object, "synonyms"))
                .chromosome(JsonFunctions.nullableString(object, "chromosome"))
                .mapLocation(JsonFunctions.nullableString(object, "mapLocation"))
                .descriptions(extractDescriptions(object.getAsJsonArray("geneDescriptions")))
                .canonicalTranscript(JsonFunctions.nullableString(object, "canonicalTranscript"))
                .geneRole(JsonFunctions.string(object, "geneRole"))
                .createDate(DateConverter.toDate(JsonFunctions.string(object, "createDate")))
                .updateDate(DateConverter.toDate(JsonFunctions.nullableString(object, "updateDate")))
                .clinicalTrials(extractClinicalTrials(object.getAsJsonArray("clinicalTrials")))
                .evidence(extractEvidence(object.getAsJsonArray("evidence")))
                .variants(extractVariants(object.getAsJsonArray("variants")))
                .molecularProfiles(extractMolecularProfiles(object.getAsJsonArray("molecularProfiles")))
                .categoryVariants(extractCategoryVariants(object.getAsJsonArray("categoryVariants")))
                .build();
    }

    @NotNull
    private static List<DescriptionInfo> extractDescriptions(@NotNull JsonArray jsonArray) {
        List<DescriptionInfo> descriptions = Lists.newArrayList();
        JsonDatamodelChecker descriptionChecker = GeneDataModelChecker.descriptionObjectChecker();

        for (JsonElement description : jsonArray) {
            JsonObject descriptionJsonObject = description.getAsJsonObject();
            descriptionChecker.check(descriptionJsonObject);

            descriptions.add(ImmutableDescriptionInfo.builder()
                    .description(JsonFunctions.string(descriptionJsonObject, "description"))
                    .references(extractReferences(descriptionJsonObject.getAsJsonArray("references")))
                    .build());
        }
        return descriptions;
    }

    @NotNull
    private static List<ReferenceInfo> extractReferences(@NotNull JsonArray jsonArray) {
        List<ReferenceInfo> references = Lists.newArrayList();
        JsonDatamodelChecker referenceChecker = GeneDataModelChecker.referenceObjectChecker();

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
    private static List<ClinicalTrialInfo> extractClinicalTrials(@NotNull JsonArray jsonArray) {
        List<ClinicalTrialInfo> clinicalTrials = Lists.newArrayList();
        JsonDatamodelChecker clinicalTrialChecker = GeneDataModelChecker.clinicalTrialObjectChecker();

        for (JsonElement clinicalTrial : jsonArray) {
            JsonObject clinicalTrialJsonObject = clinicalTrial.getAsJsonObject();
            clinicalTrialChecker.check(clinicalTrialJsonObject);

            if (JsonFunctions.nullableString(clinicalTrialJsonObject, "phase") == null) {
                LOGGER.warn("phase of study '{}' is nullable from GeneReader", JsonFunctions.string(clinicalTrialJsonObject, "nctId"));
            }
            clinicalTrials.add(ImmutableClinicalTrialInfo.builder()
                    .nctId(JsonFunctions.string(clinicalTrialJsonObject, "nctId"))
                    .title(JsonFunctions.string(clinicalTrialJsonObject, "title"))
                    .phase(JsonFunctions.nullableString(clinicalTrialJsonObject, "phase"))
                    .recruitment(JsonFunctions.string(clinicalTrialJsonObject, "recruitment"))
                    .therapies(extractTherapies(clinicalTrialJsonObject.getAsJsonArray("therapies")))
                    .build());
        }
        return clinicalTrials;
    }

    @NotNull
    private static List<TherapyInfo> extractTherapies(@NotNull JsonArray jsonArray) {
        List<TherapyInfo> therapies = Lists.newArrayList();
        JsonDatamodelChecker therapyChecker = GeneDataModelChecker.therapyObjectChecker();

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
    private static List<EvidenceInfo> extractEvidence(@NotNull JsonArray jsonArray) {
        List<EvidenceInfo> evidences = Lists.newArrayList();
        JsonDatamodelChecker evidenceChecker = GeneDataModelChecker.evidenceObjectChecker();

        for (JsonElement evidence : jsonArray) {
            JsonObject evidenceObject = evidence.getAsJsonObject();
            evidenceChecker.check(evidenceObject);

            evidences.add(ImmutableEvidenceInfo.builder()
                    .id(JsonFunctions.integer(evidenceObject, "id"))
                    .approvalStatus(JsonFunctions.string(evidenceObject, "approvalStatus"))
                    .evidenceType(JsonFunctions.string(evidenceObject, "evidenceType"))
                    .efficacyEvidence(JsonFunctions.string(evidenceObject, "efficacyEvidence"))
                    .molecularProfile(extractMolecularProfileObject(evidenceObject.getAsJsonObject("molecularProfile")))
                    .therapy(extractTherapyObject(evidenceObject.getAsJsonObject("therapy")))
                    .indication(extractIndicationObject(evidenceObject.getAsJsonObject("indication")))
                    .responseType(JsonFunctions.string(evidenceObject, "responseType"))
                    .references(extractReferences(evidenceObject.getAsJsonArray("references")))
                    .ampCapAscoEvidenceLevel(JsonFunctions.string(evidenceObject, "ampCapAscoEvidenceLevel"))
                    .ampCapAscoInferredTier(JsonFunctions.string(evidenceObject, "ampCapAscoInferredTier"))
                    .build());
        }
        return evidences;
    }

    @NotNull
    private static MolecularProfileInfo extractMolecularProfileObject(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker molecularProfileChecker = GeneDataModelChecker.molecularProfileObjectChecker();
        molecularProfileChecker.check(jsonObject);

        return ImmutableMolecularProfileInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .profileName(JsonFunctions.string(jsonObject, "profileName"))
                .build();
    }

    @NotNull
    private static TherapyInfo extractTherapyObject(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker therapyChecker = GeneDataModelChecker.therapyObjectChecker();
        therapyChecker.check(jsonObject);

        return ImmutableTherapyInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .therapyName(JsonFunctions.string(jsonObject, "therapyName"))
                .synonyms(JsonFunctions.optionalStringList(jsonObject, "synonyms"))
                .build();
    }

    @NotNull
    private static IndicationInfo extractIndicationObject(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker indicationChecker = GeneDataModelChecker.indicationObjectChecker();
        indicationChecker.check(jsonObject);

        return ImmutableIndicationInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .name(JsonFunctions.string(jsonObject, "name"))
                .source(JsonFunctions.string(jsonObject, "source"))
                .build();
    }

    @NotNull
    private static List<VariantInfo> extractVariants(@NotNull JsonArray jsonArray) {
        List<VariantInfo> variants = Lists.newArrayList();
        JsonDatamodelChecker variantChecker = GeneDataModelChecker.variantObjectChecker();

        for (JsonElement variant : jsonArray) {
            JsonObject variantObject = variant.getAsJsonObject();
            variantChecker.check(variantObject);

            variants.add(ImmutableVariantInfo.builder()
                    .id(JsonFunctions.integer(variantObject, "id"))
                    .fullName(JsonFunctions.string(variantObject, "fullName"))
                    .impact(JsonFunctions.nullableString(variantObject, "impact"))
                    .proteinEffect(JsonFunctions.nullableString(variantObject, "proteinEffect"))
                    .descriptions(extractVariantDescriptions(variantObject.getAsJsonArray("geneVariantDescriptions")))
                    .build());
        }
        return variants;
    }

    @NotNull
    private static List<DescriptionInfo> extractVariantDescriptions(@NotNull JsonArray jsonArray) {
        List<DescriptionInfo> variantDescriptions = Lists.newArrayList();
        JsonDatamodelChecker variantDescriptionChecker = GeneDataModelChecker.variantDescriptionObjectChecker();

        for (JsonElement variantDescription : jsonArray) {
            JsonObject variantDescriptionObject = variantDescription.getAsJsonObject();
            variantDescriptionChecker.check(variantDescriptionObject);

            variantDescriptions.add(ImmutableDescriptionInfo.builder()
                    .description(JsonFunctions.string(variantDescriptionObject, "description"))
                    .references(extractReferences(variantDescriptionObject.getAsJsonArray("references")))
                    .build());
        }
        return variantDescriptions;
    }

    @NotNull
    private static List<MolecularProfileInfo> extractMolecularProfiles(@NotNull JsonArray jsonArray) {
        List<MolecularProfileInfo> molecularProfiles = Lists.newArrayList();
        JsonDatamodelChecker molecularProfileChecker = GeneDataModelChecker.molecularProfileObjectChecker();

        for (JsonElement molecularProfile : jsonArray) {
            JsonObject molecularProfileObject = molecularProfile.getAsJsonObject();
            molecularProfileChecker.check(molecularProfileObject);

            molecularProfiles.add(ImmutableMolecularProfileInfo.builder()
                    .id(JsonFunctions.integer(molecularProfileObject, "id"))
                    .profileName(JsonFunctions.string(molecularProfileObject, "profileName"))
                    .treatmentApproaches(extractProfileTreatmentApproaches(molecularProfileObject.getAsJsonArray(
                            "profileTreatmentApproaches")))
                    .build());
        }
        return molecularProfiles;
    }

    @NotNull
    private static List<TreatmentApproachInfo> extractProfileTreatmentApproaches(@NotNull JsonArray jsonArray) {
        List<TreatmentApproachInfo> profileTreatmentApproaches = Lists.newArrayList();
        JsonDatamodelChecker profileTreatmentApproachChecker = GeneDataModelChecker.profileTreatmentApproachObjectChecker();

        for (JsonElement profileTreatmentApproach : jsonArray) {
            JsonObject profileTreatmentApproachObject = profileTreatmentApproach.getAsJsonObject();
            profileTreatmentApproachChecker.check(profileTreatmentApproachObject);

            profileTreatmentApproaches.add(ImmutableTreatmentApproachInfo.builder()
                    .id(JsonFunctions.integer(profileTreatmentApproachObject, "id"))
                    .name(JsonFunctions.string(profileTreatmentApproachObject, "name"))
                    .profileName(JsonFunctions.string(profileTreatmentApproachObject, "profileName"))
                    .build());
        }
        return profileTreatmentApproaches;
    }

    @NotNull
    private static List<VariantInfo> extractCategoryVariants(@NotNull JsonArray jsonArray) {
        List<VariantInfo> categoryVariants = Lists.newArrayList();
        JsonDatamodelChecker categoryVariantObjectChecker = GeneDataModelChecker.categoryVariantObjectChecker();

        for (JsonElement categoryVariant : jsonArray) {
            JsonObject categoryVariantObject = categoryVariant.getAsJsonObject();
            categoryVariantObjectChecker.check(categoryVariantObject);

            categoryVariants.add(ImmutableVariantInfo.builder()
                    .id(JsonFunctions.integer(categoryVariantObject, "id"))
                    .fullName(JsonFunctions.string(categoryVariantObject, "fullName"))
                    .impact(JsonFunctions.nullableString(categoryVariantObject, "impact"))
                    .proteinEffect(JsonFunctions.nullableString(categoryVariantObject, "proteinEffect"))
                    .descriptions(extractVariantDescriptions(categoryVariantObject.getAsJsonArray("geneVariantDescriptions")))
                    .build());
        }
        return categoryVariants;
    }
}
