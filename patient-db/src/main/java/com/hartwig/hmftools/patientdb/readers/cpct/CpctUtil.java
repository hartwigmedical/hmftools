package com.hartwig.hmftools.patientdb.readers.cpct;

import java.util.Map;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ecrf.EcrfModel;
import com.hartwig.hmftools.common.ecrf.reader.CodeList;
import com.hartwig.hmftools.common.ecrf.reader.Item;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class CpctUtil {

    private static final Logger LOGGER = LogManager.getLogger(CpctUtil.class);

    private static final String FIELD_HOSPITAL1 = "FLD.ELIGIBILITY.HOSPITAL";
    private static final String FIELD_HOSPITAL2 = "FLD.SELCRIT.NHOSPITAL";

    private CpctUtil() {
    }

    @NotNull
    public static Map<Integer, String> extractHospitalMap(@NotNull final EcrfModel model) {
        final Map<Integer, String> hospitals = Maps.newHashMap();

        final Map<String, CodeList> codeLists = model.datamodel().codeLists();

        final Item hospitalItem1 = model.datamodel().items().get(FIELD_HOSPITAL1);
        if (hospitalItem1 != null) {
            hospitals.putAll(codeLists.get(hospitalItem1.codeListOID()).values());
        } else {
            LOGGER.warn("Could not find hospital item in datamodel: " + FIELD_HOSPITAL1);
        }

        final Item hospitalItem2 = model.datamodel().items().get(FIELD_HOSPITAL2);
        if (hospitalItem2 != null) {
            hospitals.putAll(codeLists.get(hospitalItem2.codeListOID()).values());
        } else {
            LOGGER.warn("Could not find hospital item in datamodel: " + FIELD_HOSPITAL2);
        }

        return ImmutableMap.copyOf(hospitals);
    }
}
