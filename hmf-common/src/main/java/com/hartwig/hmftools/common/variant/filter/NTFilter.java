package com.hartwig.hmftools.common.variant.filter;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;

public class NTFilter implements VariantContextFilter
{
    @Override
    public boolean test(final VariantContext record)
    {
        // If NT attribute is present we want to only include ref calls (no het/hom/conflict).
        return !record.getCommonInfo().hasAttribute("NT") || record.getCommonInfo().getAttribute("NT").equals("ref");
    }
}
