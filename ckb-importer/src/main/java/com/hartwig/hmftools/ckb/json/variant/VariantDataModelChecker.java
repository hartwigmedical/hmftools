package com.hartwig.hmftools.ckb.json.variant;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;

import org.jetbrains.annotations.NotNull;

final class VariantDataModelChecker {

    private VariantDataModelChecker(){
    }

    @NotNull
    public static JsonDatamodelChecker variantObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("fullName", true);
        map.put("impact", true);
        map.put("proteinEffect", true);
        map.put("geneVariantDescriptions", true);
        map.put("type", true);
        map.put("gene", true);
        map.put("variant", true);
        map.put("createDate", true);
        map.put("updateDate", true);
        map.put("referenceTranscriptCoordinates", true);
        map.put("partnerGenes", true);
        map.put("categoryVariantPaths", true);
        map.put("evidence", true);
        map.put("extendedEvidence", true);
        map.put("molecularProfiles", true);
        map.put("allTranscriptCoordinates", true);
        map.put("memberVariants", true);
        map.put("associatedWithDrugResistance", false); //TODO implement
        map.put("transformingActivity", false); //TODO implement
        map.put("polymorphism", false); //TODO implement

        return new JsonDatamodelChecker("VariantObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker geneVariantDescriptionObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("description", true);
        map.put("references", true);

        return new JsonDatamodelChecker("VariantGeneVariantDescriptionObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker referenceObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("pubMedId", true);
        map.put("title", true);
        map.put("url", true);

        return new JsonDatamodelChecker("VariantReferenceObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker geneObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("geneSymbol", true);
        map.put("terms", true);

        return new JsonDatamodelChecker("VariantGeneObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker referenceTranscriptCoordinateObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("transcript", true);
        map.put("gDna", true);
        map.put("cDna", true);
        map.put("protein", true);
        map.put("sourceDb", true);
        map.put("refGenomeBuild", true);

        return new JsonDatamodelChecker("VariantReferenceTranscriptCoordinateObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker partnerGeneObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("gene", true);

        return new JsonDatamodelChecker("VariantPartnerGeneObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker categoryVariantPathObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("variantPath", true);
        map.put("variants", true);

        return new JsonDatamodelChecker("VariantCategoryVariantPathObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker variantVariantObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("fullName", true);
        map.put("impact", true);
        map.put("proteinEffect", true);
        map.put("associatedWithDrugResistance", false); //TODO implement
        map.put("transformingActivity", false); //TODO implement
        map.put("polymorphism", false); //TODO implement

        return new JsonDatamodelChecker("VariantVariantObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker evidenceObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("approvalStatus", true);
        map.put("evidenceType", true);
        map.put("efficacyEvidence", true);
        map.put("molecularProfile", true);
        map.put("therapy", true);
        map.put("indication", true);
        map.put("responseType", true);
        map.put("references", true);
        map.put("ampCapAscoEvidenceLevel", true);
        map.put("ampCapAscoInferredTier", true);

        return new JsonDatamodelChecker("VariantEvidenceObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker molecularProfileObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("profileName", true);
        map.put("profileTreatmentApproaches", false);

        return new JsonDatamodelChecker("VariantMolecularProfileObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker therapyChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("therapyName", true);
        map.put("synonyms", true);

        return new JsonDatamodelChecker("VariantTherapyObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker indicationChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("name", true);
        map.put("source", true);

        return new JsonDatamodelChecker("VariantIndicationObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker extendedEvidenceObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("approvalStatus", true);
        map.put("evidenceType", true);
        map.put("efficacyEvidence", true);
        map.put("molecularProfile", true);
        map.put("therapy", true);
        map.put("indication", true);
        map.put("responseType", true);
        map.put("references", true);
        map.put("ampCapAscoEvidenceLevel", true);
        map.put("ampCapAscoInferredTier", true);

        return new JsonDatamodelChecker("VariantExtendedEvidenceObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker profileTreatmentApproachObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("name", true);
        map.put("profileName", true);

        return new JsonDatamodelChecker("VariantProfileTreatmentApproachObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker allTranscriptCoordinateObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("transcript", true);
        map.put("gDna", true);
        map.put("cDna", true);
        map.put("protein", true);
        map.put("sourceDb", true);
        map.put("refGenomeBuild", true);

        return new JsonDatamodelChecker("VariantAllTranscriptCoordinateObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker memberVariantObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("fullName", true);
        map.put("impact", true);
        map.put("proteinEffect", true);
        map.put("geneVariantDescriptions", true);

        return new JsonDatamodelChecker("VariantMemberVariantObject", map);
    }
}
