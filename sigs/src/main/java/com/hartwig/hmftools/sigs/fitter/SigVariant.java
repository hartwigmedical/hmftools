package com.hartwig.hmftools.sigs.fitter;

import com.hartwig.hmftools.common.variant.SmallVariant;

import htsjdk.variant.variantcontext.VariantContext;

public class SigVariant
{
    public final SmallVariant Variant;
    public final VariantContext VariantContext;

    public SigVariant(final SmallVariant variant, final VariantContext variantContext)
    {
        Variant = variant;
        VariantContext = variantContext;
    }

    public String toString() { return Variant.toString(); }
}
