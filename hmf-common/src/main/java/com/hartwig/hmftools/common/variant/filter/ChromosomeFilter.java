package com.hartwig.hmftools.common.variant.filter;

import com.hartwig.hmftools.common.chromosome.HumanChromosome;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;

public class ChromosomeFilter implements VariantContextFilter {

    @Override
    public boolean test(final VariantContext record) {
        return HumanChromosome.contains(record.getContig());
    }
}
