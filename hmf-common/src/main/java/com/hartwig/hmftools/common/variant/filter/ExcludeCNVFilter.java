package com.hartwig.hmftools.common.variant.filter;

import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;

public class ExcludeCNVFilter implements VariantContextFilter {
    @Override
    public boolean test(final VariantContext variantContext) {
        return variantContext.getStructuralVariantType() == null || !variantContext.getStructuralVariantType()
                .equals(StructuralVariantType.CNV);
    }
}
