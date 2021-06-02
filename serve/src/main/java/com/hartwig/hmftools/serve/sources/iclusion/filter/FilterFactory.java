package com.hartwig.hmftools.serve.sources.iclusion.filter;

import java.util.Set;

import com.google.common.collect.Sets;

final class FilterFactory {

    static final Set<String> MUTATIONS_TO_FILTER = Sets.newHashSet();

    static final Set<FilterKey> MUTATION_KEYS_TO_FILTER = Sets.newHashSet();

    static {
        populateMutationsToFilter();
        populateMutationKeysToFilter();
    }

    private static void populateMutationsToFilter() {
        // We don't interpret "expression" evidence on WGS
        MUTATIONS_TO_FILTER.add("EXPRESSION");

        // We don't consider wild-type events yet.
        MUTATIONS_TO_FILTER.add("Wild-type");
    }

    private static void populateMutationKeysToFilter() {
        // Exon 12 for NPM1 not present in canonical transcript
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("NPM1", "EXON 12 MUTATION"));

        // Not a valid mutation string
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("ERBB2", "P780INS"));

        // Mutations which are inconsistent with our current gene panel
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("FLT1", "AMPLIFICATION"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("MST1R", "AMPLIFICATION"));

        MUTATION_KEYS_TO_FILTER.add(new FilterKey("FLT1", "OVEREXPRESSION"));

        MUTATION_KEYS_TO_FILTER.add(new FilterKey("ERBB4", "ACTIVATING MUTATION"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("FLT1", "ACTIVATING MUTATION"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("MAP2K4", "ACTIVATING MUTATION"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("MAP3K1", "ACTIVATING MUTATION"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("MST1R", "ACTIVATING MUTATION"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("NOTCH1", "ACTIVATING MUTATION"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("NOTCH2", "ACTIVATING MUTATION"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("NOTCH3", "ACTIVATING MUTATION"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("NOTCH4", "ACTIVATING MUTATION"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("NRG1", "ACTIVATING MUTATION"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("PIK3R1", "ACTIVATING MUTATION"));

        // Fusions that would not get reported anyways.
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("FGFR4", "FUSION"));
    }

    private FilterFactory() {
    }
}

