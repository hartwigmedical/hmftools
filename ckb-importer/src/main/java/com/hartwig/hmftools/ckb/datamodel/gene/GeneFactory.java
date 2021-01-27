package com.hartwig.hmftools.ckb.datamodel.gene;

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
import com.hartwig.hmftools.ckb.datamodel.common.ClinicalTrialInfo;
import com.hartwig.hmftools.ckb.datamodel.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.datamodel.common.EffectInfo;
import com.hartwig.hmftools.ckb.datamodel.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableClinicalTrialInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableDescriptionInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableEffectInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableEvidenceInfo;
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
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;
import com.hartwig.hmftools.common.utils.json.JsonFunctions;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class GeneFactory {

    private static final Logger LOGGER = LogManager.getLogger(GeneFactory.class);

    private GeneFactory() {

    }

    @NotNull
    public static List<Gene> readingGenes(@NotNull String geneDir) throws IOException {
        LOGGER.info("Start reading genes");

        List<Gene> genes = Lists.newArrayList();
        File[] filesGenes = new File(geneDir).listFiles();

        if (filesGenes != null) {
            LOGGER.info("The total files in the genes dir is {}", filesGenes.length);

            for (File gene : filesGenes) {
                JsonParser parser = new JsonParser();
                JsonReader reader = new JsonReader(new FileReader(gene));
                reader.setLenient(true);

                while (reader.peek() != JsonToken.END_DOCUMENT) {
                    JsonObject geneEntryObject = parser.parse(reader).getAsJsonObject();
                    JsonDatamodelChecker geneChecker = GeneDataModelChecker.geneObjectChecker();
                    geneChecker.check(geneEntryObject);

                    genes.add(ImmutableGene.builder()
                            .id(JsonFunctions.string(geneEntryObject, "id"))
                            .geneSymbol(JsonFunctions.string(geneEntryObject, "geneSymbol"))
                            .term(JsonFunctions.stringList(geneEntryObject, "terms"))
                            .entrezId(JsonFunctions.nullableString(geneEntryObject, "entrezId"))
                            .synonym(JsonFunctions.stringList(geneEntryObject, "synonyms"))
                            .chromosome(JsonFunctions.nullableString(geneEntryObject, "chromosome"))
                            .mapLocation(JsonFunctions.nullableString(geneEntryObject, "mapLocation"))
                            .description(extractGeneDescriptions(geneEntryObject.getAsJsonArray("geneDescriptions")))
                            .canonicalTranscript(JsonFunctions.nullableString(geneEntryObject, "canonicalTranscript"))
                            .geneRole(JsonFunctions.string(geneEntryObject, "geneRole"))
                            .createDate(JsonFunctions.string(geneEntryObject, "createDate"))
                            .updateDate(JsonFunctions.nullableString(geneEntryObject, "updateDate"))
                            .clinicalTrial(extractGeneClinicalTrial(geneEntryObject.getAsJsonArray("clinicalTrials")))
                            .evidence(extractGeneEvidence(geneEntryObject.getAsJsonArray("evidence")))
                            .variant(extractGeneVariant(geneEntryObject.getAsJsonArray("variants")))
                            .molecularProfile(extarctMolecularProfile(geneEntryObject.getAsJsonArray("molecularProfiles")))
                            .categoryVariant(extractCategoryVariant(geneEntryObject.getAsJsonArray("categoryVariants")))
                            .build());
                }
                reader.close();
            }
        }
        LOGGER.info("Finished reading genes");

        return genes;
    }

    @NotNull
    public static List<DescriptionInfo> extractGeneDescriptions(@NotNull JsonArray jsonArray) {
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
    public static List<ReferenceInfo> extractGeneReferences(@NotNull JsonArray jsonArray) {
        List<ReferenceInfo> references = Lists.newArrayList();
        JsonDatamodelChecker geneReferenceChecker = GeneDataModelChecker.geneReferenceObjectChecker();

        for (JsonElement geneReference : jsonArray) {
            JsonObject geneReferenceJsonObject = geneReference.getAsJsonObject();
            geneReferenceChecker.check(geneReferenceJsonObject);

            references.add(ImmutableReferenceInfo.builder()
                    .id(JsonFunctions.string(geneReferenceJsonObject, "id"))
                    .pubMedId(JsonFunctions.nullableString(geneReferenceJsonObject, "pubMedId"))
                    .title(JsonFunctions.nullableString(geneReferenceJsonObject, "title"))
                    .url(JsonFunctions.nullableString(geneReferenceJsonObject, "url"))
                    .build());

        }
        return references;
    }

    @NotNull
    public static List<ClinicalTrialInfo> extractGeneClinicalTrial(@NotNull JsonArray jsonArray) {
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
    public static List<TherapyInfo> extractGeneTherapy(@NotNull JsonArray jsonArray) {
        List<TherapyInfo> geneTherapies = Lists.newArrayList();
        JsonDatamodelChecker geneTherapiesChecker = GeneDataModelChecker.geneTherapiesObjectChecker();

        for (JsonElement geneTherapy : jsonArray) {
            JsonObject geneTherapyJsonObject = geneTherapy.getAsJsonObject();
            geneTherapiesChecker.check(geneTherapyJsonObject);

            geneTherapies.add(ImmutableTherapyInfo.builder()
                    .id(JsonFunctions.string(geneTherapyJsonObject, "id"))
                    .therapyName(JsonFunctions.string(geneTherapyJsonObject, "therapyName"))
                    .synonyms(JsonFunctions.nullableString(geneTherapyJsonObject, "synonyms"))
                    .build());
        }
        return geneTherapies;
    }

    @NotNull
    public static List<EvidenceInfo> extractGeneEvidence(@NotNull JsonArray jsonArray) {
        List<EvidenceInfo> geneEvidences = Lists.newArrayList();
        JsonDatamodelChecker geneEvidenceChecker = GeneDataModelChecker.geneEvidenceObjectChecker();

        for (JsonElement geneEvidence : jsonArray) {
            JsonObject geneEvidenceObject = geneEvidence.getAsJsonObject();
            geneEvidenceChecker.check(geneEvidenceObject);

            geneEvidences.add(ImmutableEvidenceInfo.builder()
                    .id(JsonFunctions.string(geneEvidenceObject, "id"))
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
    public static MolecularProfileInfo extractGeneMolecularProfileObject(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker geneMolecularProfileChecker = GeneDataModelChecker.geneMolecularProfileObjectChecker();
        geneMolecularProfileChecker.check(jsonObject);

        return ImmutableMolecularProfileInfo.builder()
                .id(JsonFunctions.string(jsonObject, "id"))
                .profileName(JsonFunctions.string(jsonObject, "profileName"))
                .build();
    }

    @NotNull
    public static TherapyInfo extractGeneTherapyObject(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker geneTherapyChecker = GeneDataModelChecker.geneTherapyObjectChecker();
        geneTherapyChecker.check(jsonObject);

        return ImmutableTherapyInfo.builder()
                .id(JsonFunctions.string(jsonObject, "id"))
                .therapyName(JsonFunctions.string(jsonObject, "therapyName"))
                .synonyms(JsonFunctions.nullableString(jsonObject, "synonyms"))
                .build();
    }

    @NotNull
    public static IndicationInfo extractGeneIndicationObject(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker geneIndicationChecker = GeneDataModelChecker.geneIndicationObjectChecker();
        geneIndicationChecker.check(jsonObject);

        return ImmutableIndicationInfo.builder()
                .id(JsonFunctions.string(jsonObject, "id"))
                .name(JsonFunctions.string(jsonObject, "name"))
                .source(JsonFunctions.string(jsonObject, "source"))
                .build();
    }

    @NotNull
    public static List<VariantInfo> extractGeneVariant(@NotNull JsonArray jsonArray) {
        List<VariantInfo> geneVariants = Lists.newArrayList();
        JsonDatamodelChecker geneVariantChecker = GeneDataModelChecker.geneVariantObjectChecker();

        for (JsonElement geneVariant : jsonArray) {
            JsonObject geneVariantObject = geneVariant.getAsJsonObject();
            geneVariantChecker.check(geneVariantObject);

            geneVariants.add(ImmutableVariantInfo.builder()
                    .id(JsonFunctions.string(geneVariantObject, "id"))
                    .fullName(JsonFunctions.string(geneVariantObject, "fullName"))
                    .impact(JsonFunctions.nullableString(geneVariantObject, "impact"))
                    .proteinEffect(JsonFunctions.nullableString(geneVariantObject, "proteinEffect"))
                    .description(extractVariantDescription(geneVariantObject.getAsJsonArray("geneVariantDescriptions")))
                    .build());
        }
        return geneVariants;
    }

    @NotNull
    public static List<DescriptionInfo> extractVariantDescription(@NotNull JsonArray jsonArray) {
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
    public static List<MolecularProfileInfo> extarctMolecularProfile(@NotNull JsonArray jsonArray) {
        List<MolecularProfileInfo> molecularProfiles = Lists.newArrayList();
        JsonDatamodelChecker geneMolecularProfileChecker = GeneDataModelChecker.geneMolecularProfileObjectChecker();

        for (JsonElement molecularProfile : jsonArray) {
            JsonObject molecularProfileObject = molecularProfile.getAsJsonObject();
            geneMolecularProfileChecker.check(molecularProfileObject);

            molecularProfiles.add(ImmutableMolecularProfileInfo.builder()
                    .id(JsonFunctions.string(molecularProfileObject, "id"))
                    .profileName(JsonFunctions.string(molecularProfileObject, "profileName"))
                    .treatmentApproach(extractProfileTreatmentApproach(molecularProfileObject.getAsJsonArray("profileTreatmentApproaches")))
                    .build());
        }
        return molecularProfiles;
    }

    @NotNull
    public static List<TreatmentApproachInfo> extractProfileTreatmentApproach(@NotNull JsonArray jsonArray) {
        List<TreatmentApproachInfo> geneProfileTreatmentApproaches = Lists.newArrayList();
        JsonDatamodelChecker geneProfileTreatmentApprochChecker = GeneDataModelChecker.geneProfileTreatmentApproachObjectChecker();

        for (JsonElement profileTreatmentApproach : jsonArray) {
            JsonObject profileTreatmentApproachObject = profileTreatmentApproach.getAsJsonObject();
            geneProfileTreatmentApprochChecker.check(profileTreatmentApproachObject);

            geneProfileTreatmentApproaches.add(ImmutableTreatmentApproachInfo.builder()
                    .id(JsonFunctions.string(profileTreatmentApproachObject, "id"))
                    .name(JsonFunctions.string(profileTreatmentApproachObject, "name"))
                    .profileName(JsonFunctions.string(profileTreatmentApproachObject, "profileName"))
                    .build());
        }
        return geneProfileTreatmentApproaches;
    }

    @NotNull
    public static List<EffectInfo> extractCategoryVariant(@NotNull JsonArray jsonArray) {
        List<EffectInfo> geneCategoryVariants = Lists.newArrayList();
        JsonDatamodelChecker geneProfileTreatmentApprochChecker = GeneDataModelChecker.geneCategoryVariantObjectChecker();

        for (JsonElement geneCategoryVariant : jsonArray) {
            JsonObject geneCategoryVariantObject = geneCategoryVariant.getAsJsonObject();
            geneProfileTreatmentApprochChecker.check(geneCategoryVariantObject);

            geneCategoryVariants.add(ImmutableEffectInfo.builder()
                    .id(JsonFunctions.string(geneCategoryVariantObject, "id"))
                    .fullName(JsonFunctions.string(geneCategoryVariantObject, "fullName"))
                    .impact(JsonFunctions.nullableString(geneCategoryVariantObject, "impact"))
                    .proteinEffect(JsonFunctions.nullableString(geneCategoryVariantObject, "proteinEffect"))
                    .description(extractVariantDescription(geneCategoryVariantObject.getAsJsonArray("geneVariantDescriptions")))
                    .build());
        }
        return geneCategoryVariants;
    }


}
