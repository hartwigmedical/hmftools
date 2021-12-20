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

    }

    private static void populateMutationKeysToFilter() {
        // Exon 12 for NPM1 not present in canonical transcript
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("NPM1", "EXON 12 MUTATION"));

        // Events on genes while genes not present in driver catalog
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("MRE11", "INACTIVATING MUTATION"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("MRE11", "LOSS"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("ABL1", "BCR-ABL"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("ABL1", "BCR-ABL T315I"));

        // Not a valid mutation string
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("ERBB2", "P780INS"));

        // Mutations which are inconsistent with our current gene panel
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("FLT1", "AMPLIFICATION"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("MST1R", "AMPLIFICATION"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("ERBB4", "COPY-GAIN"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("CSF1", "OVEREXPRESSION"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("FLT1", "OVEREXPRESSION"));

        MUTATION_KEYS_TO_FILTER.add(new FilterKey("ERBB4", "ACTIVATING MUTATION"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("FLT1", "ACTIVATING MUTATION"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("MAP2K4", "ACTIVATING MUTATION"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("MAP3K1", "ACTIVATING MUTATION"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("MST1R", "ACTIVATING MUTATION"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("JAK1", "ACTIVATING MUTATION"));

        // Fusions that would not get reported anyways.
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("FGFR4", "FUSION"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("RSPO4", "FUSION"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("RSPO1", "FUSION"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("ERBB4", "FUSION"));
        MUTATION_KEYS_TO_FILTER.add(new FilterKey("ERBB2", "FUSION"));
    }

    private FilterFactory() {
    }
}