package com.hartwig.hmftools.ckb.json.reference;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;

import org.jetbrains.annotations.NotNull;

final class ReferenceDataModelChecker {

    private ReferenceDataModelChecker(){
    }

    @NotNull
    public static JsonDatamodelChecker referenceObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("pubMedId", true);
        map.put("title", true);
        map.put("url", true);
        map.put("authors", true);
        map.put("journal", true);
        map.put("volume", true);
        map.put("issue", true);
        map.put("date", true);
        map.put("abstractText", true);
        map.put("year", true);
        map.put("drugs", true);
        map.put("genes", true);
        map.put("evidence", true);
        map.put("therapies", true);
        map.put("treatmentApproaches", true);
        map.put("variants", true);

        return new JsonDatamodelChecker("ReferenceObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker drugObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("drugName", true);
        map.put("terms", true);

        return new JsonDatamodelChecker("ReferenceDrugObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker geneObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("geneSymbol", true);
        map.put("terms", true);

        return new JsonDatamodelChecker("ReferenceGeneObject", map);
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

        return new JsonDatamodelChecker("ReferenceEvidenceObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker molecularProfileObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("profileName", true);

        return new JsonDatamodelChecker("ReferenceMolecularProfileObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker therapyChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("therapyName", true);
        map.put("synonyms", true);

        return new JsonDatamodelChecker("ReferenceTherapyObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker indicationChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("name", true);
        map.put("source", true);

        return new JsonDatamodelChecker("ReferenceIndicationObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker referenceChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("pubMedId", true);
        map.put("title", true);
        map.put("url", true);

        return new JsonDatamodelChecker("ReferenceReferenceObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker treatmentApproachObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("name", true);
        map.put("profileName", true);

        return new JsonDatamodelChecker("ReferenceTreatmentApproachObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker variantObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("fullName", true);
        map.put("impact", true);
        map.put("proteinEffect", true);

        return new JsonDatamodelChecker("ReferenceVariantObject", map);
    }
}
