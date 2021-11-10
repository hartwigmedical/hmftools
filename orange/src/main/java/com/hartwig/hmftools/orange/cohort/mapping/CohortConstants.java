package com.hartwig.hmftools.orange.cohort.mapping;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

public final class CohortConstants {

    public static final String PAN_CANCER_COHORT = "Pan-cancer";
    public static final String OTHER_COHORT = "Other";

    public static final List<Set<String>> DOID_COMBINATIONS_TO_MAP_TO_OTHER = Lists.newArrayList();

    static {
        // The combination of liver cancer and bile duct is occasionally used but not easily map.
        DOID_COMBINATIONS_TO_MAP_TO_OTHER.add(Sets.newHashSet("686", "4947"));
    }

    private CohortConstants() {
    }
}
