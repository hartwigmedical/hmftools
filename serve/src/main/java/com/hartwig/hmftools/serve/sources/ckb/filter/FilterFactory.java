package com.hartwig.hmftools.serve.sources.ckb.filter;

import java.util.Set;

import com.google.common.collect.Sets;

final class FilterFactory {

    static final Set<String> VARIANT_KEYWORDS_TO_FILTER = Sets.newHashSet();

    static final Set<String> GENES_FOR_WHICH_TO_FILTER_ALL = Sets.newHashSet();
    static final Set<String> GENES_FOR_WHICH_TO_FILTER_EXON_EVENTS = Sets.newHashSet();

    static {
        populateVariantKeywordsToFilter();
        populateGenesForWhichToFilterAll();
        populateGenesForWhichToFilterExonEvents();
    }

    private static void populateVariantKeywordsToFilter() {
        // We don't consider wild-type events yet.
        VARIANT_KEYWORDS_TO_FILTER.add("wild-type");

        // We cannot determine methylation with WGS/WTS
        VARIANT_KEYWORDS_TO_FILTER.add("hypermethylation");

        // This is a variant which serves as a placeholder for evidence not related to any variant.
        VARIANT_KEYWORDS_TO_FILTER.add("unknown");

        // "Expression" is not observed on DNA level
        VARIANT_KEYWORDS_TO_FILTER.add("positive");

        // We don't consider LOH a driver on its own
        VARIANT_KEYWORDS_TO_FILTER.add("LOH");
    }

    private static void populateGenesForWhichToFilterAll() {
        // COX2 lies on MT and we don't handle that in hmftools.
        GENES_FOR_WHICH_TO_FILTER_ALL.add("COX2");
    }

    private static void populateGenesForWhichToFilterExonEvents() {
        // NPM1 generally uses a different transcript, so exon events are too risky to interpret.
        GENES_FOR_WHICH_TO_FILTER_EXON_EVENTS.add("NPM1");
    }

    private FilterFactory() {
    }
}
