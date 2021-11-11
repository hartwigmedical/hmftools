package com.hartwig.hmftools.linx.germline;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PON_FILTER_PON;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;

public class GermlineFilter implements VariantContextFilter
{
    public static final String GERMLINE_MIN_QUAL = "germline_min_qual";

    private static final String FILTER_MIN_QUAL = "minQual";

    private final int mMinQualLimit;

    public GermlineFilter(int minQualLimit)
    {
        mMinQualLimit = minQualLimit;
    }

    @Override
    public boolean test(final VariantContext record)
    {
        if(record.getFilters().isEmpty())
            return true;

        if(mMinQualLimit > 0)
        {
            if(record.getPhredScaledQual() < mMinQualLimit)
                return false;

            // accept PASS, PON or 'PON;minQual' or minQual
            return record.getFilters().stream().allMatch(x -> x.equals(PASS) || x.equals(PON_FILTER_PON) || x.equals(FILTER_MIN_QUAL));
        }
        else
        {
            if(record.getFilters().size() > 1)
                return false;

            return record.getFilters().stream().anyMatch(x -> x.equals(PASS) || x.equals(PON_FILTER_PON));
        }
    }
}
