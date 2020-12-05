package com.hartwig.hmftools.linx.types;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.PON_FILTER_PON;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;

public class GermlineFilter implements VariantContextFilter
{
    private final boolean mAllowPon;

    public GermlineFilter(boolean allowPon)
    {
        mAllowPon = allowPon;
    }

    @Override
    public boolean test(final VariantContext record)
    {
        if(record.getFilters().isEmpty())
            return true;

        if(record.getFilters().size() > 1)
            return false;

        if(record.getFilters().stream().anyMatch(x -> x.equals(PASS)))
            return true;

        if(mAllowPon && record.getFilters().stream().anyMatch(x -> x.equals(PON_FILTER_PON)))
            return true;

        return false;
    }
}
