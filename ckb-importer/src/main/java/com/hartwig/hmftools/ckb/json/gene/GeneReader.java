package com.hartwig.hmftools.ckb.json.gene;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.ckb.json.CkbJsonDirectoryReader;
import com.hartwig.hmftools.ckb.json.common.ClinicalTrialInfo;
import com.hartwig.hmftools.ckb.json.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.json.common.EffectInfo;
import com.hartwig.hmftools.ckb.json.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableClinicalTrialInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableDescriptionInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableEffectInfo;
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

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class GeneReader extends CkbJsonDirectoryReader<Gene> {

    public GeneReader(@Nullable final Integer maxFilesToRead) {
        super(maxFilesToRead);
    }

    @NotNull
    @Override
    protected Gene read(@NotNull final JsonObject object) {
        JsonDatamodelChecker geneChecker = GeneDataModelChecker.geneObjectChecker();
        geneChecker.check(object);

        return ImmutableGene.builder()
                .id(JsonFunctions.integer(object, "id"))
                .geneSymbol(JsonFunctions.string(object, "geneSymbol"))
                .term(JsonFunctions.stringList(object, "terms"))
                .entrezId(JsonFunctions.nullableString(object, "entrezId"))
                .synonym(JsonFunctions.stringList(object, "synonyms"))
                .chromosome(JsonFunctions.nullableString(object, "chromosome"))
                .mapLocation(JsonFunctions.nullableString(object, "mapLocation"))
                .description(extractGeneDescriptions(object.getAsJsonArray("geneDescriptions")))
                .canonicalTranscript(JsonFunctions.nullableString(object, "canonicalTranscript"))
                .geneRole(JsonFunctions.string(object, "geneRole"))
                .createDate(DateConverter.toDate(JsonFunctions.string(object, "createDate")))
                .updateDate(DateConverter.toDate(JsonFunctions.nullableString(object, "updateDate")))
                .clinicalTrial(extractGeneClinicalTrial(object.getAsJsonArray("clinicalTrials")))
                .evidence(extractGeneEvidence(object.getAsJsonArray("evidence")))
                .variant(extractGeneVariant(object.getAsJsonArray("variants")))
                .molecularProfile(extarctMolecularProfile(object.getAsJsonArray("molecularProfiles")))
                .categoryVariant(extractCategoryVariant(object.getAsJsonArray("categoryVariants")))
                .build();
    }

    @NotNull
    private static List<DescriptionInfo> extractGeneDescriptions(@NotNull JsonArray jsonArray) {
        List<DescriptionInfo> geneDescriptions = Lists.newArrayList();
        JsonDatamodelChecker geneDescriptionChecker = GeneDataModelChecker.geneDescriptionObjectChecker();

        for (JsonElement geneDescription : jsonArray) {
            JsonObject geneDescriptionJsonObject = geneDescription.getAsJsonObject();
            geneDescriptionChecker.check(geneDescriptionJsonObject);

            geneDescriptions.add(ImmutableDescriptionInfo.builder()
                    .description(JsonFunctions.string(geneDescriptionJsonObject, "description"))
                    .reference(extractGeneReferences(geneDescriptionJsonObject.getAsJsonArray("references")))
                    .build());
        }
        return geneDescriptions;
    }

    @NotNull
    private static List<ReferenceInfo> extractGeneReferences(@NotNull JsonArray jsonArray) {
        List<ReferenceInfo> references = Lists.newArrayList();
        JsonDatamodelChecker geneReferenceChecker = GeneDataModelChecker.geneReferenceObjectChecker();

        for (JsonElement geneReference : jsonArray) {
            JsonObject geneReferenceJsonObject = geneReference.getAsJsonObject();
            geneReferenceChecker.check(geneReferenceJsonObject);

            references.add(ImmutableReferenceInfo.builder()
                    .id(JsonFunctions.integer(geneReferenceJsonObject, "id"))
                    .pubMedId(JsonFunctions.nullableString(geneReferenceJsonObject, "pubMedId"))
                    .title(JsonFunctions.nullableString(geneReferenceJsonObject, "title"))
                    .url(JsonFunctions.nullableString(geneReferenceJsonObject, "url"))
                    .build());

        }
        return references;
    }

    @NotNull
    private static List<ClinicalTrialInfo> extractGeneClinicalTrial(@NotNull JsonArray jsonArray) {
        List<ClinicalTrialInfo> clinicalTrials = Lists.newArrayList();
        JsonDatamodelChecker geneClinicalTrialChecker = GeneDataModelChecker.geneClinicalTrialObjectChecker();

        for (JsonElement geneClinicalTrial : jsonArray) {
            JsonObject geneClinicalTrialJsonObject = geneClinicalTrial.getAsJsonObject();
            geneClinicalTrialChecker.check(geneClinicalTrialJsonObject);

            clinicalTrials.add(ImmutableClinicalTrialInfo.builder()
                    .nctId(JsonFunctions.string(geneClinicalTrialJsonObject, "nctId"))
                    .title(JsonFunctions.string(geneClinicalTrialJsonObject, "title"))
                    .phase(JsonFunctions.string(geneClinicalTrialJsonObject, "phase"))
                    .recruitment(JsonFunctions.string(geneClinicalTrialJsonObject, "recruitment"))
                    .therapy(extractGeneTherapy(geneClinicalTrialJsonObject.getAsJsonArray("therapies")))
                    .build());
        }
        return clinicalTrials;
    }

    @NotNull
    private static List<TherapyInfo> extractGeneTherapy(@NotNull JsonArray jsonArray) {
        List<TherapyInfo> geneTherapies = Lists.newArrayList();
        JsonDatamodelChecker geneTherapiesChecker = GeneDataModelChecker.geneTherapiesObjectChecker();

        for (JsonElement geneTherapy : jsonArray) {
            JsonObject geneTherapyJsonObject = geneTherapy.getAsJsonObject();
            geneTherapiesChecker.check(geneTherapyJsonObject);

            geneTherapies.add(ImmutableTherapyInfo.builder()
                    .id(JsonFunctions.integer(geneTherapyJsonObject, "id"))
                    .therapyName(JsonFunctions.string(geneTherapyJsonObject, "therapyName"))
                    .synonyms(JsonFunctions.nullableString(geneTherapyJsonObject, "synonyms"))
                    .build());
        }
        return geneTherapies;
    }

    @NotNull
    private static List<EvidenceInfo> extractGeneEvidence(@NotNull JsonArray jsonArray) {
        List<EvidenceInfo> geneEvidences = Lists.newArrayList();
        JsonDatamodelChecker geneEvidenceChecker = GeneDataModelChecker.geneEvidenceObjectChecker();

        for (JsonElement geneEvidence : jsonArray) {
            JsonObject geneEvidenceObject = geneEvidence.getAsJsonObject();
            geneEvidenceChecker.check(geneEvidenceObject);

            geneEvidences.add(ImmutableEvidenceInfo.builder()
                    .id(JsonFunctions.integer(geneEvidenceObject, "id"))
                    .approvalStatus(JsonFunctions.string(geneEvidenceObject, "approvalStatus"))
                    .evidenceType(JsonFunctions.string(geneEvidenceObject, "evidenceType"))
                    .efficacyEvidence(JsonFunctions.string(geneEvidenceObject, "efficacyEvidence"))
                    .molecularProfile(extractGeneMolecularProfileObject(geneEvidenceObject.getAsJsonObject("molecularProfile")))
                    .therapy(extractGeneTherapyObject(geneEvidenceObject.getAsJsonObject("therapy")))
                    .indication(extractGeneIndicationObject(geneEvidenceObject.getAsJsonObject("indication")))
                    .responseType(JsonFunctions.string(geneEvidenceObject, "responseType"))
                    .reference(extractGeneReferences(geneEvidenceObject.getAsJsonArray("references")))
                    .ampCapAscoEvidenceLevel(JsonFunctions.string(geneEvidenceObject, "ampCapAscoEvidenceLevel"))
                    .ampCapAscoInferredTier(JsonFunctions.string(geneEvidenceObject, "ampCapAscoInferredTier"))
                    .build());
        }
        return geneEvidences;
    }

    @NotNull
    private static MolecularProfileInfo extractGeneMolecularProfileObject(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker geneMolecularProfileChecker = GeneDataModelChecker.geneMolecularProfileObjectChecker();
        geneMolecularProfileChecker.check(jsonObject);

        return ImmutableMolecularProfileInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .profileName(JsonFunctions.string(jsonObject, "profileName"))
                .build();
    }

    @NotNull
    private static TherapyInfo extractGeneTherapyObject(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker geneTherapyChecker = GeneDataModelChecker.geneTherapyObjectChecker();
        geneTherapyChecker.check(jsonObject);

        return ImmutableTherapyInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .therapyName(JsonFunctions.string(jsonObject, "therapyName"))
                .synonyms(JsonFunctions.nullableString(jsonObject, "synonyms"))
                .build();
    }

    @NotNull
    private static IndicationInfo extractGeneIndicationObject(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker geneIndicationChecker = GeneDataModelChecker.geneIndicationObjectChecker();
        geneIndicationChecker.check(jsonObject);

        return ImmutableIndicationInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .name(JsonFunctions.string(jsonObject, "name"))
                .source(JsonFunctions.string(jsonObject, "source"))
                .build();
    }

    @NotNull
    private static List<VariantInfo> extractGeneVariant(@NotNull JsonArray jsonArray) {
        List<VariantInfo> geneVariants = Lists.newArrayList();
        JsonDatamodelChecker geneVariantChecker = GeneDataModelChecker.geneVariantObjectChecker();

        for (JsonElement geneVariant : jsonArray) {
            JsonObject geneVariantObject = geneVariant.getAsJsonObject();
            geneVariantChecker.check(geneVariantObject);

            geneVariants.add(ImmutableVariantInfo.builder()
                    .id(JsonFunctions.integer(geneVariantObject, "id"))
                    .fullName(JsonFunctions.string(geneVariantObject, "fullName"))
                    .impact(JsonFunctions.nullableString(geneVariantObject, "impact"))
                    .proteinEffect(JsonFunctions.nullableString(geneVariantObject, "proteinEffect"))
                    .description(extractVariantDescription(geneVariantObject.getAsJsonArray("geneVariantDescriptions")))
                    .build());
        }
        return geneVariants;
    }

    @NotNull
    private static List<DescriptionInfo> extractVariantDescription(@NotNull JsonArray jsonArray) {
        List<DescriptionInfo> geneVariantDescriptions = Lists.newArrayList();
        JsonDatamodelChecker geneVariantDescriptionChecker = GeneDataModelChecker.geneVariantDescriptionObjectChecker();

        for (JsonElement variantDescription : jsonArray) {
            JsonObject variantDescriptionObject = variantDescription.getAsJsonObject();
            geneVariantDescriptionChecker.check(variantDescriptionObject);

            geneVariantDescriptions.add(ImmutableDescriptionInfo.builder()
                    .description(JsonFunctions.string(variantDescriptionObject, "description"))
                    .reference(extractGeneReferences(variantDescriptionObject.getAsJsonArray("references")))
                    .build());
        }
        return geneVariantDescriptions;
    }

    @NotNull
    private static List<MolecularProfileInfo> extarctMolecularProfile(@NotNull JsonArray jsonArray) {
        List<MolecularProfileInfo> molecularProfiles = Lists.newArrayList();
        JsonDatamodelChecker geneMolecularProfileChecker = GeneDataModelChecker.geneMolecularProfileObjectChecker();

        for (JsonElement molecularProfile : jsonArray) {
            JsonObject molecularProfileObject = molecularProfile.getAsJsonObject();
            geneMolecularProfileChecker.check(molecularProfileObject);

            molecularProfiles.add(ImmutableMolecularProfileInfo.builder()
                    .id(JsonFunctions.integer(molecularProfileObject, "id"))
                    .profileName(JsonFunctions.string(molecularProfileObject, "profileName"))
                    .treatmentApproach(extractProfileTreatmentApproach(molecularProfileObject.getAsJsonArray("profileTreatmentApproaches")))
                    .build());
        }
        return molecularProfiles;
    }

    @NotNull
    private static List<TreatmentApproachInfo> extractProfileTreatmentApproach(@NotNull JsonArray jsonArray) {
        List<TreatmentApproachInfo> geneProfileTreatmentApproaches = Lists.newArrayList();
        JsonDatamodelChecker geneProfileTreatmentApprochChecker = GeneDataModelChecker.geneProfileTreatmentApproachObjectChecker();

        for (JsonElement profileTreatmentApproach : jsonArray) {
            JsonObject profileTreatmentApproachObject = profileTreatmentApproach.getAsJsonObject();
            geneProfileTreatmentApprochChecker.check(profileTreatmentApproachObject);

            geneProfileTreatmentApproaches.add(ImmutableTreatmentApproachInfo.builder()
                    .id(JsonFunctions.integer(profileTreatmentApproachObject, "id"))
                    .name(JsonFunctions.string(profileTreatmentApproachObject, "name"))
                    .profileName(JsonFunctions.string(profileTreatmentApproachObject, "profileName"))
                    .build());
        }
        return geneProfileTreatmentApproaches;
    }

    @NotNull
    private static List<EffectInfo> extractCategoryVariant(@NotNull JsonArray jsonArray) {
        List<EffectInfo> geneCategoryVariants = Lists.newArrayList();
        JsonDatamodelChecker geneProfileTreatmentApprochChecker = GeneDataModelChecker.geneCategoryVariantObjectChecker();

        for (JsonElement geneCategoryVariant : jsonArray) {
            JsonObject geneCategoryVariantObject = geneCategoryVariant.getAsJsonObject();
            geneProfileTreatmentApprochChecker.check(geneCategoryVariantObject);

            geneCategoryVariants.add(ImmutableEffectInfo.builder()
                    .id(JsonFunctions.integer(geneCategoryVariantObject, "id"))
                    .fullName(JsonFunctions.string(geneCategoryVariantObject, "fullName"))
                    .impact(JsonFunctions.nullableString(geneCategoryVariantObject, "impact"))
                    .proteinEffect(JsonFunctions.nullableString(geneCategoryVariantObject, "proteinEffect"))
                    .description(extractVariantDescription(geneCategoryVariantObject.getAsJsonArray("geneVariantDescriptions")))
                    .build());
        }
        return geneCategoryVariants;
    }
}
