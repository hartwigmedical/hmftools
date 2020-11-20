package com.hartwig.hmftools.serve.sources.iclusion.filter;

import java.util.Set;

import com.google.common.collect.Sets;

final class FilterFactory {

    static final Set<String> MUTATIONS_TO_FILTER = Sets.newHashSet();

    static final Set<FilterKey> MUTATION_KEYS_TO_FILTER = Sets.newHashSet();

    static {
        populateMutationToFilter();
        populateMutationKeysToFilter();
    }

    private static void populateMutationToFilter() {
        // We don't interpret "expression" evidence on WGS
        MUTATIONS_TO_FILTER.add("EXPRESSION");

        // We don't consider wild-type events yet.
        MUTATIONS_TO_FILTER.add("Wild-type");
    }

    private static void populateMutationKeysToFilter() {
        // Exon 12 for NPM1 not present in canonical transcript
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("NPM1", "EXON 12 MUTATION"));

        // Not a valid mutation string
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("ERBB2 (HER2)", "P780INS"));
    }

    private FilterFactory() {
    }
}

