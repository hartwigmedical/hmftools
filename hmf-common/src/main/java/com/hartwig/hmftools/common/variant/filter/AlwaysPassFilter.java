package com.hartwig.hmftools.common.variant.filter;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;

public class AlwaysPassFilter implements VariantContextFilter {
    @Override
    public boolean test(final VariantContext variantContext) {
        return true;
    }
}
