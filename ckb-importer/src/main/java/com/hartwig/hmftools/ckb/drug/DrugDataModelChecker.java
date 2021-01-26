package com.hartwig.hmftools.ckb.drug;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;

import org.jetbrains.annotations.NotNull;

public final class DrugDataModelChecker {

    private DrugDataModelChecker() {
    }

    @NotNull
    public static JsonDatamodelChecker drugObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("drugName", true);
        map.put("terms", true);
        map.put("synonyms", true);
        map.put("tradeName", true);
        map.put("drugDescriptions", true);
        map.put("drugClasses", true);
        map.put("casRegistryNum", true);
        map.put("ncitId", true);
        map.put("createDate", true);
        map.put("clinicalTrials", true);
        map.put("evidence", true);
        map.put("therapies", true);
        map.put("globalApprovalStatus", true);

        return new JsonDatamodelChecker("DrugObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker drugDescriptionObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("description", true);
        map.put("references", true);

        return new JsonDatamodelChecker("DrugDescriptionObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker drugReferenceObjectChecker() {
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

        return new JsonDatamodelChecker("DrugReferenceObject", map);

    }

    @NotNull
    public static JsonDatamodelChecker drugClassObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("drugClass", true);

        return new JsonDatamodelChecker("DrugClassObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker drugClinicalTrialObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("nctId", true);
        map.put("title", true);
        map.put("phase", true);
        map.put("recruitment", true);
        map.put("therapies", true);

        return new JsonDatamodelChecker("DrugClinicalTrialObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker drugClinicalTrialTherapyObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("therapyName", true);
        map.put("synonyms", true);

        return new JsonDatamodelChecker("DrugClinicalTrialTherapyObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker drugEvidenceObjectChecker() {
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

        return new JsonDatamodelChecker("DrugClinicalTrialTherapyObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker drugEvidenceMolecularProfileObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("profileName", true);
        map.put("profileTreatmentApproache", false); //check if needed

        return new JsonDatamodelChecker("DrugEvidenceMolecularProfileObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker drugEvidenceTherapyObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("therapyName", true);
        map.put("synonyms", true);

        return new JsonDatamodelChecker("DrugEvidenceTherapyObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker drugEvidenceIndicationObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("name", true);
        map.put("source", true);

        return new JsonDatamodelChecker("DrugEvidenceIndicationObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker drugEvidenceReferenceObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("pubMedId", true);
        map.put("title", true);
        map.put("url", true);

        return new JsonDatamodelChecker("DrugEvidenceReferenceObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker drugTherapyObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("therapyName", true);
        map.put("synonyms", true);

        return new JsonDatamodelChecker("DrugTherapyObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker drugGlobalApprovalStatusObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("therapy", true);
        map.put("indication", true);
        map.put("molecularProfile", true);
        map.put("approvalAuthority", true);
        map.put("approvalStatus", true);

        return new JsonDatamodelChecker("DrugGlobalApprovalStatusObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker drugGlobalApprovalStatusTherapyObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("therapyName", true);
        map.put("synonyms", true);

        return new JsonDatamodelChecker("DrugGlobalApprovalStatusTherapyObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker drugGlobalApprovalStatusIndicationObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("name", true);
        map.put("source", true);

        return new JsonDatamodelChecker("DrugGlobalApprovalStatusIndicationObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker drugGlobalApprovalStatusMolecularProfileObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("profileName", false);
        map.put("profileTreatmentApproache", null); //check if needed

        return new JsonDatamodelChecker("DrugGlobalApprovalStatusMolecularProfileObject", map);
    }

}
