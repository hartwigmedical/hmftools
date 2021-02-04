package com.hartwig.hmftools.ckb.reader.indication;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;

import org.jetbrains.annotations.NotNull;

public class IndicationDataModelChecker {

    private IndicationDataModelChecker() {

    }

    @NotNull
    public static JsonDatamodelChecker indicationObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("name", true);
        map.put("source", true);
        map.put("definition", true);
        map.put("currentPreferredTerm", true);
        map.put("lastUpdateDateFromDO", true);
        map.put("altIds", true);
        map.put("termId", true);
        map.put("evidence", true);
        map.put("clinicalTrials", true);

        return new JsonDatamodelChecker("IndicationObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker indicationEvidenceObjectChecker() {
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

        return new JsonDatamodelChecker("IndicationEvidenceObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker indicationEvidenceMolecularProfileObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("profileName", true);

        return new JsonDatamodelChecker("IndicationEvidenceMolecularProfileObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker indicationEvidenceTherapyObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("therapyName", true);
        map.put("synonyms", true);

        return new JsonDatamodelChecker("IndicationEvidenceTherapyObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker indicationEvidenceIndicationObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("name", true);
        map.put("source", true);

        return new JsonDatamodelChecker("IndicationEvidenceIndicationObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker indicationEvidenceReferenceObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("pubMedId", true);
        map.put("title", true);
        map.put("url", true);

        return new JsonDatamodelChecker("IndicationEvidenceReferenceObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker indicationClinicaltrialObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("nctId", true);
        map.put("title", true);
        map.put("phase", true);
        map.put("recruitment", true);
        map.put("therapies", true);

        return new JsonDatamodelChecker("IndicationClinicalTrialObject", map);
    }
}
