package com.hartwig.hmftools.common.variant.filter;

import java.util.Set;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;

public class PassingFilter implements VariantContextFilter {

    @Override
    public boolean test(final VariantContext record) {
        Set<String> filters = record.getFilters();
        return filters.isEmpty() || (filters.size() == 1 && filters.contains("PASS"));
    }
}
