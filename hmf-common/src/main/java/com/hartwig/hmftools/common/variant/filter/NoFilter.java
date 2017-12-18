package com.hartwig.hmftools.common.variant.filter;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;

public class NoFilter implements VariantContextFilter {

    @Override
    public boolean test(final VariantContext record) {
        return true;
    }
}
