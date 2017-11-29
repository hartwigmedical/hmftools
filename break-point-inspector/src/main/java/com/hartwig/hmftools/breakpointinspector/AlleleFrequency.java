package com.hartwig.hmftools.breakpointinspector;

import org.apache.commons.lang3.tuple.Pair;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

class AlleleFrequency {

    static final String VCF_INFO_TAG = "BPI_AF";

    static void UpdateVCFHeader(final VCFHeader header) {
        header.addMetaDataLine(new VCFInfoHeaderLine(VCF_INFO_TAG, 2, VCFHeaderLineType.Float, "AF at each breakpoint"));
    }

    private static Double calculate(final BreakpointStats bp, final double support) {
        final double total = bp.PR_Only_Normal + bp.PR_SR_Normal + support;
        return total > 0.0 ? support / total : 0.0;
    }

    static Pair<Double, Double> calculate(final SampleStats stats) {
        final BreakpointStats bp1 = stats.BP1_Stats, bp2 = stats.BP2_Stats;
        final double support = bp1.PR_SR_Support + bp1.PR_Only_Support + bp1.SR_Only_Support + bp2.SR_Only_Support;
        return Pair.of(calculate(stats.BP1_Stats, support), calculate(stats.BP2_Stats, support));
    }

}
