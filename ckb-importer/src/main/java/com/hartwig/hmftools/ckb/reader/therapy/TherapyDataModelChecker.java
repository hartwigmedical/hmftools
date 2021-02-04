package com.hartwig.hmftools.ckb.reader.therapy;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;

import org.jetbrains.annotations.NotNull;

public class TherapyDataModelChecker {

    private TherapyDataModelChecker() {

    }

    @NotNull
    public static JsonDatamodelChecker therapyObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("therapyName", true);
        map.put("synonyms", true);
        map.put("therapyDescriptions", true);
        map.put("createDate", true);
        map.put("updateDate", true);
        map.put("evidence", true);
        map.put("clinicalTrials", true);
        map.put("drugs", true);
        map.put("globalApprovalStatus", true);

        return new JsonDatamodelChecker("TherapyObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker descriptionObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("description", true);
        map.put("references", true);

        return new JsonDatamodelChecker("TherapyDescriptionObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker referencesObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("pubMedId", true);
        map.put("title", true);
        map.put("url", true);

        return new JsonDatamodelChecker("TherapyReferencesObject", map);
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

        return new JsonDatamodelChecker("TherapyEvidenceObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker molecularProfileObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("profileName", true);

        return new JsonDatamodelChecker("TherapyMolecularprofileObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker therapyChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("therapyName", true);
        map.put("synonyms", true);

        return new JsonDatamodelChecker("TherapyTherapyObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker indicationChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("name", true);
        map.put("source", true);

        return new JsonDatamodelChecker("TherapyIndicationObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker referenceChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("pubMedId", true);
        map.put("title", true);
        map.put("url", true);

        return new JsonDatamodelChecker("TherapyReferenceObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker clinicalTrialChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("nctId", true);
        map.put("title", true);
        map.put("phase", true);
        map.put("recruitment", true);
        map.put("therapies", true);

        return new JsonDatamodelChecker("TherapyClinicalTrialObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker drugsChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("drugName", true);
        map.put("terms", true);

        return new JsonDatamodelChecker("TherapyDrugslObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker globalApprovalStatusChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("therapy", true);
        map.put("indication", true);
        map.put("molecularProfile", true);
        map.put("approvalAuthority", true);
        map.put("approvalStatus", true);

        return new JsonDatamodelChecker("TherapyGlobalApprovalStatusObject", map);
    }
}
