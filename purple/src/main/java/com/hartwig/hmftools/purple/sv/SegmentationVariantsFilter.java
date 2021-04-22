package com.hartwig.hmftools.purple.sv;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.INFERRED;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;

public class SegmentationVariantsFilter extends PassingVariantFilter {

    @Override
    public boolean test(final VariantContext variantContext) {
        return super.test(variantContext) || variantContext.hasAttribute(INFERRED);
    }
}
