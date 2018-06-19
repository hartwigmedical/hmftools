package com.hartwig.hmftools.breakpointinspector;

import com.hartwig.hmftools.breakpointinspector.datamodel.BreakpointStats;
import com.hartwig.hmftools.breakpointinspector.datamodel.SampleStats;

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
        final double total = bp.prOnlyNormal + bp.prSrNormal + support;
        return total > 0.0 ? support / total : 0.0;
    }

    @NotNull
    static Pair<Double, Double> calculate(@NotNull final SampleStats stats) {
        final BreakpointStats bp1 = stats.bp1Stats, bp2 = stats.bp2Stats;
        final double support = bp1.prSrSupport + bp1.prOnlySupport + bp1.srOnlySupport + bp2.srOnlySupport;
        return Pair.of(calculate(stats.bp1Stats, support), calculate(stats.bp2Stats, support));
    }
}
