package com.hartwig.hmftools.ckb.reader.molecularprofile;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;

import org.jetbrains.annotations.NotNull;

public class MolecularProfileDataModelChecker {

    private MolecularProfileDataModelChecker() {

    }

    @NotNull
    public static JsonDatamodelChecker molecularProfileObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("profileName", true);
        map.put("geneVariants", true);
        map.put("profileTreatmentApproaches", true);
        map.put("createDate", true);
        map.put("updateDate", true);
        map.put("complexMolecularProfileEvidence", true);
        map.put("treatmentApproachEvidence", true);
        map.put("variantAssociatedClinicalTrials", true);
        map.put("variantLevelEvidence", true);
        map.put("extendedEvidence", true);

        return new JsonDatamodelChecker("MolecularProfileObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker molecularProfileGeneVariantObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("fullName", true);
        map.put("impact", true);
        map.put("proteinEffect", true);

        return new JsonDatamodelChecker("MolecularProfileGeneVariantObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker molecularProfileProfileTreatmentApproacjObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("name", true);
        map.put("profileName", true);

        return new JsonDatamodelChecker("MolecularProfileProfileTreatmentApproachObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker molecularProfileComplexMolecularProfileEvidence() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("totalCount", true);
        map.put("complexMolecularProfileEvidence", true);

        return new JsonDatamodelChecker("MolecularProfileComplexMolecularProfileEvidenceObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker molecularProfileComplexMolecularProfileEvidenceList() {
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
        map.put("relevantTreatmentApproaches", true);

        return new JsonDatamodelChecker("MolecularProfileComplexMolecularProfileEvidenceListObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker molecularProfileMolecularprofile() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("profileName", true);

        return new JsonDatamodelChecker("MolecularProfileMolecularprofileObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker molecularProfileTherapy() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("therapyName", true);
        map.put("synonyms", true);

        return new JsonDatamodelChecker("MolecularProfileTherapyObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker molecularProfileIndication() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("name", true);
        map.put("source", true);

        return new JsonDatamodelChecker("MolecularProfileIndicationObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker molecularProfileReference() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("pubMedId", true);
        map.put("title", true);
        map.put("url", true);

        return new JsonDatamodelChecker("MolecularProfileReferenceObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker molecularProfileRelevantTreatmentApproach() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("name", true);
        map.put("profileName", true);

        return new JsonDatamodelChecker("MolecularProfileRelevantTreatmentApproachObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker molecularProfileTreatmentApproachEvidence() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("totalCount", true);
        map.put("treatmentApproachEvidence", true);

        return new JsonDatamodelChecker("MolecularProfileTreatmentApproachEvidenceObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker molecularProfileTreatmentApproachEvidenceList() {
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
        map.put("relevantTreatmentApproaches", true);

        return new JsonDatamodelChecker("MolecularProfileTreatmentApproachEvidenceListObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker variantAssociatedClinicalTrial() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("nctId", true);
        map.put("title", true);
        map.put("phase", true);
        map.put("recruitment", true);
        map.put("therapies", true);

        return new JsonDatamodelChecker("MolecularProfileVariantAssociatedClinicalTrialObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker molecularProfilevariantLevelEvidence() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("totalCount", true);
        map.put("variantLevelEvidences", true);

        return new JsonDatamodelChecker("MolecularProfileVariantLevelEvidenceObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker molecularProfilevariantLevelEvidenceList() {
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
        map.put("relevantTreatmentApproaches", true);

        return new JsonDatamodelChecker("MolecularProfileVariantLevelEvidenceObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker extendedEvidenceEvidence() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("totalCount", true);
        map.put("extendedEvidence", true);

        return new JsonDatamodelChecker("MolecularProfileExtendedEvidenceObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker extendedEvidenceList() {
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

        return new JsonDatamodelChecker("MolecularProfileExtendedEvidenceListObject", map);
    }
}
