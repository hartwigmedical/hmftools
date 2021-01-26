package com.hartwig.hmftools.ckb.globaltherapyapprovalstatus;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;

import org.jetbrains.annotations.NotNull;

public class GlobalTherapyApprovalStatusDataModelChecker {

    private GlobalTherapyApprovalStatusDataModelChecker() {

    }

    @NotNull
    public static JsonDatamodelChecker globalTherapyApprovalStatusObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("totalCount", true);
        map.put("globalTherapyApprovalStatuses", true);

        return new JsonDatamodelChecker("GlobalTherapyApprovalStatusObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker globalTherapyApprovalStatusListObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("therapy", true);
        map.put("indication", true);
        map.put("molecularProfile", true);
        map.put("approvalAuthority", true);
        map.put("approvalStatus", true);

        return new JsonDatamodelChecker("GlobalTherapyApprovalStatusListObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker globalTherapyApprovalStatusTherapyObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("therapyName", true);
        map.put("synonyms", true);

        return new JsonDatamodelChecker("GlobalTherapyApprovalStatusTherapyObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker globalTherapyApprovalStatusIndicationObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("name", true);
        map.put("source", true);

        return new JsonDatamodelChecker("GlobalTherapyApprovalStatusIndicationObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker globalTherapyApprovalStatusMolecularprofileObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("profileName", true);
        map.put("profileTreatmentApproache", false); //check if needed

        return new JsonDatamodelChecker("GlobalTherapyApprovalStatusMolecularProfileObject", map);
    }

}
