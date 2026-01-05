package com.hartwig.hmftools.linx.germline;

import static com.hartwig.hmftools.common.sv.SvVcfTags.PON_FILTER_PON;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS_FILTER;

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

        return record.getFilters().stream().anyMatch(x -> x.equals(PASS_FILTER) || x.equals(PON_FILTER_PON));
    }
}
