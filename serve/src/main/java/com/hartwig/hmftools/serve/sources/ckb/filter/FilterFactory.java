package com.hartwig.hmftools.serve.sources.ckb.filter;

import java.util.Set;

import com.google.common.collect.Sets;

public class FilterFactory {

    static final Set<String> PROFILE_NAME_TO_FILTER = Sets.newHashSet();

    static {
        populateMutationToFilter();
    }

    private static void populateMutationToFilter() {

        // We don't consider wild-type events yet.
        PROFILE_NAME_TO_FILTER.add("wild-type");

        // We cannot determine methylation with WGS/WTS
        PROFILE_NAME_TO_FILTER.add("hypermethylation");

        // Is not clear what profile it is
        PROFILE_NAME_TO_FILTER.add("unknown");

        // "Expression" is not observed on DNA leve
        PROFILE_NAME_TO_FILTER.add("Positive");

    }


    private FilterFactory() {
    }
}
