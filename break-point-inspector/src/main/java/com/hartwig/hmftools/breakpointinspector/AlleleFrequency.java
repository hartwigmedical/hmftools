package com.hartwig.hmftools.breakpointinspector;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

final class AlleleFrequency {

    static final String VCF_AF_INFO_TAG = "BPI_AF";

    private AlleleFrequency() {
    }

    static void updateVCFHeader(@NotNull final VCFHeader header) {
        header.addMetaDataLine(new VCFInfoHeaderLine(VCF_AF_INFO_TAG, 2, VCFHeaderLineType.Float, "AF at each breakpoint"));
    }

    @NotNull
    private static Double calculate(@NotNull final BreakpointStats bp, final double support) {
        final double total = bp.PR_Only_Normal + bp.PR_SR_Normal + support;
        return total > 0.0 ? support / total : 0.0;
    }

    @NotNull
    static Pair<Double, Double> calculate(@NotNull final SampleStats stats) {
        final BreakpointStats bp1 = stats.BP1_Stats, bp2 = stats.BP2_Stats;
        final double support = bp1.PR_SR_Support + bp1.PR_Only_Support + bp1.SR_Only_Support + bp2.SR_Only_Support;
        return Pair.of(calculate(stats.BP1_Stats, support), calculate(stats.BP2_Stats, support));
    }
}
