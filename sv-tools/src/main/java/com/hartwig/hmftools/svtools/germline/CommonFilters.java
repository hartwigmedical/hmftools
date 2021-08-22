package com.hartwig.hmftools.svtools.germline;

import static com.hartwig.hmftools.svtools.germline.FilterConstants.POLY_A;
import static com.hartwig.hmftools.svtools.germline.FilterConstants.POLY_T;
import static com.hartwig.hmftools.svtools.germline.FilterConstants.SHORT_RESCUE_LENGTH;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.HOMSEQ;

import htsjdk.variant.variantcontext.VariantContext;

public class CommonFilters
{
    public static boolean tooShortToRescue(int length) { return length < SHORT_RESCUE_LENGTH; }

    public static boolean isPolyATSequence(final VariantContext variant)
    {
        final String homology = variant.getAttributeAsString(HOMSEQ, "");
        return homology.contains(POLY_A) || homology.contains(POLY_T);
    }

}
