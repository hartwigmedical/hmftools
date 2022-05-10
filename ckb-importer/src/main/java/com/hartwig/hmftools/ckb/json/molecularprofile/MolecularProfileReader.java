package com.hartwig.hmftools.ckb.json.molecularprofile;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.ckb.json.CkbJsonDirectoryReader;
import com.hartwig.hmftools.ckb.json.common.ClinicalTrialInfo;
import com.hartwig.hmftools.ckb.json.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableClinicalTrialInfo;
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

public class MolecularProfileReader extends CkbJsonDirectoryReader<JsonMolecularProfile> {
    private static final Logger LOGGER = LogManager.getLogger(MolecularProfileReader.class);

    public MolecularProfileReader(@Nullable final Integer maxFilesToRead) {
        super(maxFilesToRead);
    }

    @NotNull
    @Override
    protected JsonMolecularProfile read(@NotNull final JsonObject object) {
        JsonDatamodelChecker molecularProfileChecker = MolecularProfileDataModelChecker.molecularProfileObjectChecker();
        molecularProfileChecker.check(object);

        return ImmutableJsonMolecularProfile.builder()
                .id(JsonFunctions.integer(object, "id"))
                .profileName(JsonFunctions.string(object, "profileName"))
                .geneVariants(extractGeneVariants(object.getAsJsonArray("geneVariants")))
                .treatmentApproaches(extractProfileTreatmentApproaches(object.getAsJsonArray("profileTreatmentApproaches")))
                .createDate(DateConverter.toDate(JsonFunctions.string(object, "createDate")))
                .updateDate(DateConverter.toDate(JsonFunctions.string(object, "updateDate")))
                .complexMolecularProfileEvidence(extractComplexMolecularProfileEvidence(object.getAsJsonObject(
                        "complexMolecularProfileEvidence")))
                .treatmentApproachEvidence(extractTreatmentApproachEvidence(object.getAsJsonObject("treatmentApproachEvidence")))
                .variantAssociatedClinicalTrials(extractVariantAssociatedClinicalTrials(object.getAsJsonArray(
                        "variantAssociatedClinicalTrials")))
                .variantLevelEvidence(extractVariantLevelEvidence(object.getAsJsonObject("variantLevelEvidence")))
                .extendedEvidence(extractExtendedEvidence(object.getAsJsonObject("extendedEvidence")))
                .build();
    }

    @NotNull
    private static List<VariantInfo> extractGeneVariants(@NotNull JsonArray jsonArray) {
        List<VariantInfo> geneVariants = Lists.newArrayList();
        JsonDatamodelChecker geneVariantChecker = MolecularProfileDataModelChecker.geneVariantObjectChecker();

        for (JsonElement geneVariant : jsonArray) {
            JsonObject geneVariantJsonObject = geneVariant.getAsJsonObject();
            geneVariantChecker.check(geneVariantJsonObject);

            geneVariants.add(ImmutableVariantInfo.builder()
                    .id(JsonFunctions.integer(geneVariantJsonObject, "id"))
                    .fullName(JsonFunctions.string(geneVariantJsonObject, "fullName"))
                    .impact(JsonFunctions.nullableString(geneVariantJsonObject, "impact"))
                    .proteinEffect(JsonFunctions.nullableString(geneVariantJsonObject, "proteinEffect"))
                    .build());
        }
        return geneVariants;
    }

    @NotNull
    private static List<TreatmentApproachInfo> extractProfileTreatmentApproaches(@NotNull JsonArray jsonArray) {
        List<TreatmentApproachInfo> profileTreatmentApproaches = Lists.newArrayList();
        JsonDatamodelChecker profileTreatmentApproachChecker = MolecularProfileDataModelChecker.profileTreatmentApproachObjectChecker();

        for (JsonElement profileTreatmentApproach : jsonArray) {
            JsonObject profileTreatmentApproachJsonObject = profileTreatmentApproach.getAsJsonObject();
            profileTreatmentApproachChecker.check(profileTreatmentApproachJsonObject);

            profileTreatmentApproaches.add(ImmutableTreatmentApproachInfo.builder()
                    .id(JsonFunctions.integer(profileTreatmentApproachJsonObject, "id"))
                    .name(JsonFunctions.string(profileTreatmentApproachJsonObject, "name"))
                    .profileName(JsonFunctions.string(profileTreatmentApproachJsonObject, "profileName"))
                    .build());
        }
        return profileTreatmentApproaches;
    }

    @NotNull
    private static JsonMolecularProfileExtendedEvidence extractComplexMolecularProfileEvidence(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker complexMolecularProfileEvidenceChecker =
                MolecularProfileDataModelChecker.complexMolecularProfileEvidenceChecker();
        complexMolecularProfileEvidenceChecker.check(jsonObject);

        return ImmutableJsonMolecularProfileExtendedEvidence.builder()
                .totalCount(JsonFunctions.integer(jsonObject, "totalCount"))
                .evidences(extractComplexMolecularProfileEvidenceList(jsonObject.getAsJsonArray("complexMolecularProfileEvidence")))
                .build();
    }

    @NotNull
    private static List<EvidenceInfo> extractComplexMolecularProfileEvidenceList(@NotNull JsonArray jsonArray) {
        List<EvidenceInfo> complexMolecularProfileEvidenceList = Lists.newArrayList();
        JsonDatamodelChecker complexMolecularProfileEvidenceListChecker =
                MolecularProfileDataModelChecker.complexMolecularProfileEvidenceListChecker();

        for (JsonElement complexMolecularProfileEvidence : jsonArray) {
            JsonObject complexMolecularProfileEvidenceJsonObject = complexMolecularProfileEvidence.getAsJsonObject();
            complexMolecularProfileEvidenceListChecker.check(complexMolecularProfileEvidenceJsonObject);

            complexMolecularProfileEvidenceList.add(ImmutableEvidenceInfo.builder()
                    .id(JsonFunctions.integer(complexMolecularProfileEvidenceJsonObject, "id"))
                    .approvalStatus(JsonFunctions.string(complexMolecularProfileEvidenceJsonObject, "approvalStatus"))
                    .evidenceType(JsonFunctions.string(complexMolecularProfileEvidenceJsonObject, "evidenceType"))
                    .efficacyEvidence(JsonFunctions.string(complexMolecularProfileEvidenceJsonObject, "efficacyEvidence"))
                    .molecularProfile(extractMolecularProfile(complexMolecularProfileEvidenceJsonObject.getAsJsonObject("molecularProfile")))
                    .therapy(extractTherapy(complexMolecularProfileEvidenceJsonObject.getAsJsonObject("therapy")))
                    .indication(extractIndication(complexMolecularProfileEvidenceJsonObject.getAsJsonObject("indication")))
                    .responseType(JsonFunctions.string(complexMolecularProfileEvidenceJsonObject, "responseType"))
                    .references(extractReferences(complexMolecularProfileEvidenceJsonObject.getAsJsonArray("references")))
                    .ampCapAscoEvidenceLevel(JsonFunctions.string(complexMolecularProfileEvidenceJsonObject, "ampCapAscoEvidenceLevel"))
                    .ampCapAscoInferredTier(JsonFunctions.string(complexMolecularProfileEvidenceJsonObject, "ampCapAscoInferredTier"))
                    .treatmentApproaches(extractRelevantTreatmentApproaches(complexMolecularProfileEvidenceJsonObject.getAsJsonArray(
                            "relevantTreatmentApproaches")))
                    .build());
        }
        return complexMolecularProfileEvidenceList;
    }

    @NotNull
    private static MolecularProfileInfo extractMolecularProfile(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker molecularProfileChecker = MolecularProfileDataModelChecker.molecularProfileChecker();
        molecularProfileChecker.check(jsonObject);

        return ImmutableMolecularProfileInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .profileName(JsonFunctions.string(jsonObject, "profileName"))
                .build();
    }

    @NotNull
    private static TherapyInfo extractTherapy(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker therapyChecker = MolecularProfileDataModelChecker.therapyChecker();
        therapyChecker.check(jsonObject);

        return ImmutableTherapyInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .therapyName(JsonFunctions.string(jsonObject, "therapyName"))
                .synonyms(JsonFunctions.optionalStringList(jsonObject, "synonyms"))
                .build();
    }

    @NotNull
    private static IndicationInfo extractIndication(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker indicationChecker = MolecularProfileDataModelChecker.indicationChecker();
        indicationChecker.check(jsonObject);

        return ImmutableIndicationInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .name(JsonFunctions.string(jsonObject, "name"))
                .source(JsonFunctions.string(jsonObject, "source"))
                .build();
    }

    @NotNull
    private static List<ReferenceInfo> extractReferences(@NotNull JsonArray jsonArray) {
        List<ReferenceInfo> references = Lists.newArrayList();
        JsonDatamodelChecker referenceChecker = MolecularProfileDataModelChecker.referenceChecker();

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
    private static List<TreatmentApproachInfo> extractRelevantTreatmentApproaches(@NotNull JsonArray jsonArray) {
        List<TreatmentApproachInfo> relevantTreatmentApproaches = Lists.newArrayList();
        JsonDatamodelChecker relevantTreatmentApproachChecker = MolecularProfileDataModelChecker.relevantTreatmentApproachChecker();

        for (JsonElement relevantTreatmentApproach : jsonArray) {
            JsonObject relevantTreatmentApproachJsonObject = relevantTreatmentApproach.getAsJsonObject();
            relevantTreatmentApproachChecker.check(relevantTreatmentApproachJsonObject);

            relevantTreatmentApproaches.add(ImmutableTreatmentApproachInfo.builder()
                    .id(JsonFunctions.integer(relevantTreatmentApproachJsonObject, "id"))
                    .name(JsonFunctions.string(relevantTreatmentApproachJsonObject, "name"))
                    .profileName(JsonFunctions.string(relevantTreatmentApproachJsonObject, "profileName"))
                    .build());
        }
        return relevantTreatmentApproaches;
    }

    @NotNull
    private static JsonMolecularProfileExtendedEvidence extractTreatmentApproachEvidence(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker treatmentApproachEvidenceChecker = MolecularProfileDataModelChecker.treatmentApproachEvidenceChecker();
        treatmentApproachEvidenceChecker.check(jsonObject);

        return ImmutableJsonMolecularProfileExtendedEvidence.builder()
                .totalCount(JsonFunctions.integer(jsonObject, "totalCount"))
                .evidences(extractTreatmentApproachEvidenceList(jsonObject.getAsJsonArray("treatmentApproachEvidence")))
                .build();
    }

    @NotNull
    private static List<EvidenceInfo> extractTreatmentApproachEvidenceList(@NotNull JsonArray jsonArray) {
        List<EvidenceInfo> treatmentApproachEvidenceList = Lists.newArrayList();
        JsonDatamodelChecker treatmentApproachEvidenceListChecker = MolecularProfileDataModelChecker.treatmentApproachEvidenceListChecker();

        for (JsonElement treatmentApproachEvidence : jsonArray) {
            JsonObject treatmentApproachEvidenceJsonObject = treatmentApproachEvidence.getAsJsonObject();
            treatmentApproachEvidenceListChecker.check(treatmentApproachEvidenceJsonObject);

            treatmentApproachEvidenceList.add(ImmutableEvidenceInfo.builder()
                    .id(JsonFunctions.integer(treatmentApproachEvidenceJsonObject, "id"))
                    .approvalStatus(JsonFunctions.string(treatmentApproachEvidenceJsonObject, "approvalStatus"))
                    .evidenceType(JsonFunctions.string(treatmentApproachEvidenceJsonObject, "evidenceType"))
                    .efficacyEvidence(JsonFunctions.string(treatmentApproachEvidenceJsonObject, "efficacyEvidence"))
                    .molecularProfile(extractMolecularProfile(treatmentApproachEvidenceJsonObject.getAsJsonObject("molecularProfile")))
                    .therapy(extractTherapy(treatmentApproachEvidenceJsonObject.getAsJsonObject("therapy")))
                    .indication(extractIndication(treatmentApproachEvidenceJsonObject.getAsJsonObject("indication")))
                    .responseType(JsonFunctions.string(treatmentApproachEvidenceJsonObject, "responseType"))
                    .references(extractReferences(treatmentApproachEvidenceJsonObject.getAsJsonArray("references")))
                    .ampCapAscoEvidenceLevel(JsonFunctions.string(treatmentApproachEvidenceJsonObject, "ampCapAscoEvidenceLevel"))
                    .ampCapAscoInferredTier(JsonFunctions.string(treatmentApproachEvidenceJsonObject, "ampCapAscoInferredTier"))
                    .treatmentApproaches(extractRelevantTreatmentApproaches(treatmentApproachEvidenceJsonObject.getAsJsonArray(
                            "relevantTreatmentApproaches")))
                    .build());
        }
        return treatmentApproachEvidenceList;
    }

    @NotNull
    private static List<ClinicalTrialInfo> extractVariantAssociatedClinicalTrials(@NotNull JsonArray jsonArray) {
        List<ClinicalTrialInfo> variantAssociatedClinicalTrials = Lists.newArrayList();
        JsonDatamodelChecker variantAssociatedClinicalTrialChecker =
                MolecularProfileDataModelChecker.variantAssociatedClinicalTrialChecker();

        for (JsonElement variantAssociatedClinicalTrial : jsonArray) {
            JsonObject variantAssociatedClinicalTrialJsonObject = variantAssociatedClinicalTrial.getAsJsonObject();
            variantAssociatedClinicalTrialChecker.check(variantAssociatedClinicalTrialJsonObject);

            String nctId = JsonFunctions.string(variantAssociatedClinicalTrialJsonObject, "nctId");
            String phase = JsonFunctions.nullableString(variantAssociatedClinicalTrialJsonObject, "phase");

            if (phase == null) {
                LOGGER.warn("phase of study '{}' is null in MolecularProfileReader", nctId);
            }

            variantAssociatedClinicalTrials.add(ImmutableClinicalTrialInfo.builder()
                    .nctId(nctId)
                    .title(JsonFunctions.string(variantAssociatedClinicalTrialJsonObject, "title"))
                    .phase(phase)
                    .recruitment(JsonFunctions.string(variantAssociatedClinicalTrialJsonObject, "recruitment"))
                    .therapies(extractTherapyList(variantAssociatedClinicalTrialJsonObject.getAsJsonArray("therapies")))
                    .build());
        }
        return variantAssociatedClinicalTrials;
    }

    @NotNull
    private static List<TherapyInfo> extractTherapyList(@NotNull JsonArray jsonArray) {
        List<TherapyInfo> therapies = Lists.newArrayList();
        JsonDatamodelChecker therapyChecker = MolecularProfileDataModelChecker.therapyChecker();

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
    private static JsonMolecularProfileExtendedEvidence extractVariantLevelEvidence(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker variantLevelEvidenceChecker = MolecularProfileDataModelChecker.variantLevelEvidenceChecker();
        variantLevelEvidenceChecker.check(jsonObject);

        return ImmutableJsonMolecularProfileExtendedEvidence.builder()
                .totalCount(JsonFunctions.integer(jsonObject, "totalCount"))
                .evidences(extractVariantLevelEvidenceList(jsonObject.getAsJsonArray("variantLevelEvidences")))
                .build();
    }

    @NotNull
    private static List<EvidenceInfo> extractVariantLevelEvidenceList(@NotNull JsonArray jsonArray) {
        List<EvidenceInfo> variantLevelEvidenceList = Lists.newArrayList();
        JsonDatamodelChecker variantLevelEvidenceListChecker = MolecularProfileDataModelChecker.variantLevelEvidenceListChecker();

        for (JsonElement variantLevelEvidence : jsonArray) {
            JsonObject variantLevelEvidenceJsonObject = variantLevelEvidence.getAsJsonObject();
            variantLevelEvidenceListChecker.check(variantLevelEvidenceJsonObject);

            variantLevelEvidenceList.add(ImmutableEvidenceInfo.builder()
                    .id(JsonFunctions.integer(variantLevelEvidenceJsonObject, "id"))
                    .approvalStatus(JsonFunctions.string(variantLevelEvidenceJsonObject, "approvalStatus"))
                    .evidenceType(JsonFunctions.string(variantLevelEvidenceJsonObject, "evidenceType"))
                    .efficacyEvidence(JsonFunctions.string(variantLevelEvidenceJsonObject, "efficacyEvidence"))
                    .molecularProfile(extractMolecularProfile(variantLevelEvidenceJsonObject.getAsJsonObject("molecularProfile")))
                    .therapy(extractTherapy(variantLevelEvidenceJsonObject.getAsJsonObject("therapy")))
                    .indication(extractIndication(variantLevelEvidenceJsonObject.getAsJsonObject("indication")))
                    .responseType(JsonFunctions.string(variantLevelEvidenceJsonObject, "responseType"))
                    .references(extractReferences(variantLevelEvidenceJsonObject.getAsJsonArray("references")))
                    .ampCapAscoEvidenceLevel(JsonFunctions.string(variantLevelEvidenceJsonObject, "ampCapAscoEvidenceLevel"))
                    .ampCapAscoInferredTier(JsonFunctions.string(variantLevelEvidenceJsonObject, "ampCapAscoInferredTier"))
                    .treatmentApproaches(extractRelevantTreatmentApproaches(variantLevelEvidenceJsonObject.getAsJsonArray(
                            "relevantTreatmentApproaches")))
                    .build());
        }
        return variantLevelEvidenceList;
    }

    @NotNull
    private static JsonMolecularProfileExtendedEvidence extractExtendedEvidence(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker extendedEvidenceChecker = MolecularProfileDataModelChecker.extendedEvidenceChecker();
        extendedEvidenceChecker.check(jsonObject);

        return ImmutableJsonMolecularProfileExtendedEvidence.builder()
                .totalCount(JsonFunctions.integer(jsonObject, "totalCount"))
                .evidences(extractExtendedEvidenceList(jsonObject.getAsJsonArray("extendedEvidence")))
                .build();
    }

    @NotNull
    private static List<EvidenceInfo> extractExtendedEvidenceList(@NotNull JsonArray jsonArray) {
        List<EvidenceInfo> extendedEvidenceList = Lists.newArrayList();
        JsonDatamodelChecker extendedEvidenceListChecker = MolecularProfileDataModelChecker.extendedEvidenceListChecker();

        for (JsonElement extendedEvidence : jsonArray) {
            JsonObject extendedEvidenceJsonObject = extendedEvidence.getAsJsonObject();
            extendedEvidenceListChecker.check(extendedEvidenceJsonObject);

            extendedEvidenceList.add(ImmutableEvidenceInfo.builder()
                    .id(JsonFunctions.integer(extendedEvidenceJsonObject, "id"))
                    .approvalStatus(JsonFunctions.string(extendedEvidenceJsonObject, "approvalStatus"))
                    .evidenceType(JsonFunctions.string(extendedEvidenceJsonObject, "evidenceType"))
                    .efficacyEvidence(JsonFunctions.string(extendedEvidenceJsonObject, "efficacyEvidence"))
                    .molecularProfile(extractMolecularProfile(extendedEvidenceJsonObject.getAsJsonObject("molecularProfile")))
                    .therapy(extractTherapy(extendedEvidenceJsonObject.getAsJsonObject("therapy")))
                    .indication(extractIndication(extendedEvidenceJsonObject.getAsJsonObject("indication")))
                    .responseType(JsonFunctions.string(extendedEvidenceJsonObject, "responseType"))
                    .references(extractReferences(extendedEvidenceJsonObject.getAsJsonArray("references")))
                    .ampCapAscoEvidenceLevel(JsonFunctions.string(extendedEvidenceJsonObject, "ampCapAscoEvidenceLevel"))
                    .ampCapAscoInferredTier(JsonFunctions.string(extendedEvidenceJsonObject, "ampCapAscoInferredTier"))

                    .build());
        }
        return extendedEvidenceList;
    }
}
