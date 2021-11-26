package com.hartwig.hmftools.gripss.filters;

import static com.hartwig.hmftools.gripss.filters.FilterConstants.POLY_A;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.POLY_T;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.SHORT_RESCUE_LENGTH;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_HOMSEQ;

import htsjdk.variant.variantcontext.VariantContext;

public class CommonFilters
{
    public static boolean tooShortToRescue(int length) { return length < SHORT_RESCUE_LENGTH; }

    public static boolean isPolyATSequence(final VariantContext variant)
    {
        final String homology = variant.getAttributeAsString(VT_HOMSEQ, "");
        return homology.contains(POLY_A) || homology.contains(POLY_T);
    }

}
