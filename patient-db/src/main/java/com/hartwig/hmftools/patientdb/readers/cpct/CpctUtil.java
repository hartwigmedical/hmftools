package com.hartwig.hmftools.patientdb.readers.cpct;

import java.util.Map;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ecrf.EcrfModel;

import org.jetbrains.annotations.NotNull;

public final class CpctUtil {

    private static final String FIELD_HOSPITAL1 = "FLD.ELIGIBILITY.HOSPITAL";
    private static final String FIELD_HOSPITAL2 = "FLD.SELCRIT.NHOSPITAL";

    private CpctUtil() {
    }

    @NotNull
    public static Map<Integer, String> extractHospitalMap(@NotNull final EcrfModel model) {
        final Map<Integer, String> hospitals = Maps.newHashMap();
        hospitals.putAll(model.datamodel().codeLists().get(model.datamodel().items().get(FIELD_HOSPITAL1).codeListOID()).values());
        hospitals.putAll(model.datamodel().codeLists().get(model.datamodel().items().get(FIELD_HOSPITAL2).codeListOID()).values());
        return ImmutableMap.copyOf(hospitals);
    }
}
