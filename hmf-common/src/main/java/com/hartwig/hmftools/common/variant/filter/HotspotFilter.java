package com.hartwig.hmftools.common.variant.filter;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;

public class HotspotFilter implements VariantContextFilter {

    private static final String HOTSPOT_TAG = "HOTSPOT";

    @Override
    public boolean test(final VariantContext record) {
        return record.hasAttribute(HOTSPOT_TAG);
    }
}
