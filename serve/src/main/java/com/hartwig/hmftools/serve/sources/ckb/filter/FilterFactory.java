package com.hartwig.hmftools.serve.sources.ckb.filter;

import java.util.Set;

import com.google.common.collect.Sets;

public final class FilterFactory {

    static final Set<String> VARIANT_KEYWORDS_TO_FILTER = Sets.newHashSet();

    static {
        populateMutationToFilter();
    }

    private static void populateMutationToFilter() {
        // We don't consider wild-type events yet.
        VARIANT_KEYWORDS_TO_FILTER.add("wild-type");

        // We cannot determine methylation with WGS/WTS
        VARIANT_KEYWORDS_TO_FILTER.add("hypermethylation");

        // Is not clear what profile it is
        VARIANT_KEYWORDS_TO_FILTER.add("unknown");

        // "Expression" is not observed on DNA level
        VARIANT_KEYWORDS_TO_FILTER.add("positive");
    }

    private FilterFactory() {
    }
}
