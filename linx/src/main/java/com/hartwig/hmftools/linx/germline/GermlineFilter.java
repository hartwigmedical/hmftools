package com.hartwig.hmftools.linx.germline;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PON_FILTER_PON;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;

public class GermlineFilter implements VariantContextFilter
{
    public GermlineFilter() {}

    @Override
    public boolean test(final VariantContext record)
    {
        if(record.getFilters().isEmpty())
            return true;

        if(record.getFilters().size() > 1)
            return false;

        return record.getFilters().stream().anyMatch(x -> x.equals(PASS) || x.equals(PON_FILTER_PON));
    }
}
