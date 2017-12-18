package com.hartwig.hmftools.common.variant.filter;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;

public class LoggingVariantContextFilter implements VariantContextFilter {

    private final VariantContextFilter filter;
    private int pass;
    private int fail;

    public LoggingVariantContextFilter(@NotNull final VariantContextFilter filter) {
        this.filter = filter;
    }

    public int getPass() {
        return pass;
    }

    public int getFail() {
        return fail;
    }

    @Override
    public boolean test(final VariantContext variantContext) {
        if (filter.test(variantContext)) {
            pass++;
            return true;
        }

        fail++;
        return false;
    }
}
