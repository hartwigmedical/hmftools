package com.hartwig.hmftools.common.variant.filter;

import java.util.List;
import java.util.Set;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public final class NearIndelPonFilter {

    private static final int DISTANCE = 10;
    private static final String SOMATIC_FLAG = "SOMATIC_PON";
    private static final String GERMLINE_FLAG = "GERMLINE_PON";

    private NearIndelPonFilter() {
    }

    public static boolean isIndelNearPon(final int index, @NotNull final List<VariantContext> contexts) {
        final VariantContext subject = contexts.get(index);
        if (!subject.isIndel() || subject.isFiltered()) {
            return false;
        }

        // JOBA: Look backwards
        for (int i = index - 1; i >= 0; i--) {
            final VariantContext query = contexts.get(i);
            int queryEnd = query.getStart() + query.getReference().length() - 1 + DISTANCE;
            if (queryEnd < subject.getStart() || !query.getContig().equals(subject.getContig())) {
                break;
            }

            if (isIndelPonFiltered(query)) {
                return true;
            }
        }

        // JOBA: Look forwards
        int subjectEnd = subject.getStart() + subject.getReference().length() - 1 + DISTANCE;
        for (int i = index + 1; i < contexts.size(); i++) {
            final VariantContext query = contexts.get(i);
            if (query.getStart() > subjectEnd || !query.getContig().equals(subject.getContig())) {
                break;
            }

            if (isIndelPonFiltered(query)) {
                return true;
            }
        }

        return false;
    }

    private static boolean isIndelPonFiltered(@NotNull final VariantContext pon) {
        final Set<String> filters = pon.getFilters();
        return pon.isIndel() && (filters.contains(SOMATIC_FLAG) || filters.contains(GERMLINE_FLAG));
    }
}
