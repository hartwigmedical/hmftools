package com.hartwig.hmftools.svtools.germline;

import htsjdk.variant.variantcontext.VariantContext;

public class SvData
{
    private final VariantContext mContext;


    public SvData(final VariantContext context)
    {
        mContext = context;

    }
}
