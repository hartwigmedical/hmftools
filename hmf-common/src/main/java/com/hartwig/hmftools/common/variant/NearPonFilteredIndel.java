package com.hartwig.hmftools.common.variant;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.utils.ArrayDeck;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public final class NearPonFilteredIndel {

    private static final int DISTANCE = 10;
    private static final String SOMATIC_FLAG = "SOMATIC_PON";
    private static final String GERMLINE_FLAG = "GERMLINE_PON";

    private NearPonFilteredIndel() {
    }

    public static boolean isPonFilteredIndel(@NotNull final VariantContext variant) {
        final Set<String> filters = variant.getFilters();
        return variant.isIndel() && (filters.contains(SOMATIC_FLAG) || filters.contains(GERMLINE_FLAG));
    }

    public static boolean isNearPonFilteredIndel(final int index, @NotNull final List<VariantContext> contexts) {
        return isNearPonFilteredIndel(index, new ArrayDeck<>(contexts));
    }

    public static boolean isNearPonFilteredIndel(final int index, @NotNull final ArrayDeck<VariantContext> contexts) {
        final VariantContext subject = contexts.get(index);
        if (!subject.isIndel() || subject.isFiltered()) {
            return false;
        }

        // Look backwards
        for (int i = index - 1; i >= 0; i--) {
            final VariantContext query = contexts.get(i);
            int queryEnd = query.getStart() + query.getReference().length() - 1 + DISTANCE;
            if (queryEnd < subject.getStart() || !query.getContig().equals(subject.getContig())) {
                break;
            }

            if (isPonFilteredIndel(query)) {
                return true;
            }
        }

        // Look forwards
        int subjectEnd = subject.getStart() + subject.getReference().length() - 1 + DISTANCE;
        for (int i = index + 1; i < contexts.size(); i++) {
            final VariantContext query = contexts.get(i);
            if (query.getStart() > subjectEnd || !query.getContig().equals(subject.getContig())) {
                break;
            }

            if (isPonFilteredIndel(query)) {
                return true;
            }
        }

        return false;
    }
}
