package com.hartwig.hmftools.common.variant.filter;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;

public class SGTFilter implements VariantContextFilter {

    @Override
    public boolean test(final VariantContext record) {
        if (record.getCommonInfo().hasAttribute("SGT")) {
            final String[] sgt = record.getCommonInfo().getAttributeAsString("SGT", "GC->AT").split("->");
            if (sgt.length == 2 && sgt[0].equals(sgt[1])) {
                return false;
            }
        }
        return true;
    }
}
