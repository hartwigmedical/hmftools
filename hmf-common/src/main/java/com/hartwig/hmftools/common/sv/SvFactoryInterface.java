package com.hartwig.hmftools.common.sv;

import java.util.List;

import com.hartwig.hmftools.common.sv.gridss.GridssSvFactory;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;

public interface SvFactoryInterface
{
    void setGenotypeOrdinals(int referenceOrdinal, int tumorOrdinal);

    void addVariantContext(final VariantContext context);

    List<StructuralVariant> results();

    List<VariantContext> unmatched();

    StructuralVariant createSV(final VariantContext first, final VariantContext second);

    StructuralVariant createSingleBreakend(final VariantContext context);

    public static SvFactoryInterface buildSvFactory(final boolean useGridssVcf, final VariantContextFilter filter)
    {
        if(useGridssVcf)
        {
            return GridssSvFactory.build(filter);
        }
        else
        {
            return StructuralVariantFactory.build(filter);
        }
    }
}
