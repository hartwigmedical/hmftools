package com.hartwig.hmftools.sage.vis;

import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;
import j2html.tags.DomContent;

public record AminoAcidElements(String geneRegionLabel, DomContent ref, DomContent alt, @Nullable VariantContext variant)
{
}
