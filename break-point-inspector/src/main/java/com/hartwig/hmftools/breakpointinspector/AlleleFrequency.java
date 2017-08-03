package com.hartwig.hmftools.breakpointinspector;

import org.apache.commons.lang3.tuple.Pair;

class AlleleFrequency {

    private static Double calculate(final HMFVariantContext ctx, final BreakpointStats bp) {
        if (ctx.isShortDelete()) {
            final double support = bp.SR_Only_Support + bp.PR_SR_Support;
            final double total = bp.PR_SR_Normal + support;
            return total > 0 ? support / total : 0.0;
        }

        final double support = bp.PR_SR_Support + bp.PR_Only_Support + bp.SR_Only_Support;
        final double total = bp.PR_SR_Normal + bp.PR_Only_Normal + support;
        return total > 0 ? support / total : 0.0;
    }

    static Pair<Double, Double> calculate(final HMFVariantContext ctx, final SampleStats tumorResult) {
        return Pair.of(calculate(ctx, tumorResult.BP1_Stats), calculate(ctx, tumorResult.BP2_Stats));
    }

}
