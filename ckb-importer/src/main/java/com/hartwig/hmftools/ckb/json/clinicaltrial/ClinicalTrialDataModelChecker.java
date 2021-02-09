package com.hartwig.hmftools.ckb.json.clinicaltrial;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;

import org.jetbrains.annotations.NotNull;

final class ClinicalTrialDataModelChecker {

    private ClinicalTrialDataModelChecker() {
    }

    @NotNull
    public static JsonDatamodelChecker clinicalTrialObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("nctId", true);
        map.put("title", true);
        map.put("phase", true);
        map.put("recruitment", true);
        map.put("therapies", true);
        map.put("ageGroups", true);
        map.put("gender", true);
        map.put("variantRequirements", true);
        map.put("sponsors", true);
        map.put("updateDate", true);
        map.put("indications", true);
        map.put("variantRequirementDetails", true);
        map.put("clinicalTrialLocations", true);

        return new JsonDatamodelChecker("ClinicalTrialObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker variantRequirementDetailObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("molecularProfile", true);
        map.put("requirementType", true);
        return new JsonDatamodelChecker("ClinicalTrialVariantRequirementDetailObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker molecularProfileObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("profileName", true);
        return new JsonDatamodelChecker("ClinicalTrialMolecularProfileObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker therapyObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("therapyName", true);
        map.put("synonyms", true);

        return new JsonDatamodelChecker("ClinicalTrialTherapyObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker indicationObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("name", true);
        map.put("source", true);
        return new JsonDatamodelChecker("ClinicalTrialIndicationObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker locationObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("nctId", true);
        map.put("facility", true);
        map.put("city", true);
        map.put("country", true);
        map.put("status", true);
        map.put("state", true);
        map.put("zip", true);
        map.put("clinicalTrialContacts", true);
        return new JsonDatamodelChecker("ClinicalTrialLocationObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker contactObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("name", true);
        map.put("email", true);
        map.put("phone", true);
        map.put("phoneExt", true);
        map.put("role", true);
        return new JsonDatamodelChecker("ClinicalTrialContactObject", map);
    }
}
