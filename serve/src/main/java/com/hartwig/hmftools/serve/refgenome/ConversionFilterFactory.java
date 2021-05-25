package com.hartwig.hmftools.serve.refgenome;

import java.util.Set;

import com.google.common.collect.Sets;

final class ConversionFilterFactory {

    static final Set<String> GENES_TO_FILTER = Sets.newHashSet();

    static {
        populateGenesToFilter();
    }

    private static void populateGenesToFilter() {
        // These genes lie on a part of the ref genome that has been flipped.
        // As a result these genes flipped strand between ref genome versions and are filtered for ref genome conversion
        GENES_TO_FILTER.add("MAGEA1");
        GENES_TO_FILTER.add("PDE4DIP");
        GENES_TO_FILTER.add("NCOA4");
    }

    private ConversionFilterFactory() {
    }
}
