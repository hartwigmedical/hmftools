package com.hartwig.hmftools.common.variant.filter;

import com.hartwig.hmftools.common.variant.enrich.HotspotEnrichment;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;

public class HotspotFilter implements VariantContextFilter {

    @Override
    public boolean test(final VariantContext record) {
        return record.hasAttribute(HotspotEnrichment.HOTSPOT_FLAG);
    }
}
