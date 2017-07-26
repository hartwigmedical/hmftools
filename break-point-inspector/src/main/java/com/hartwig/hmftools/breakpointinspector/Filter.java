package com.hartwig.hmftools.breakpointinspector;

import java.util.List;
import java.util.stream.Stream;

import com.google.common.collect.Lists;

class Filter {

    static String getFilterString(final Util.HMFVariantContext ctx, final SampleStats tumorStats, final SampleStats refStats) {

        final List<String> filters = Lists.newArrayList(ctx.Filter);

        if (ctx.Type == Util.HMFVariantType.DEL && ctx.MantaBP1.ReferenceIndex == ctx.MantaBP2.ReferenceIndex
                && (ctx.MantaBP2.Position - ctx.MantaBP1.Position) < 2000) {
            // short delete logic, must have SR support
            final int tumor_SR =
                    Stream.of(tumorStats.BP1_Stats, tumorStats.BP2_Stats).mapToInt(s -> s.PR_SR_Support + s.SR_Only_Support).sum();
            if (tumor_SR == 0) {
                filters.add("HMF_SRSupportZero");
            }

            // short delete logic, must have SR support
            final int ref_SR = Stream.of(refStats.BP1_Stats, refStats.BP2_Stats).mapToInt(s -> s.PR_SR_Support + s.SR_Only_Support).sum();
            if (ref_SR > 0) {
                filters.add("HMF_SRNormalSupport");
            }
        } else {
            // we only need to check BP1 as BP1 PR+PRSR == BP2 PR+PRSR
            if (refStats.BP1_Stats.PR_Only_Support + refStats.BP1_Stats.PR_SR_Support > 0) {
                filters.add("HMF_PRNormalSupport");
            }
        }

        boolean concordance = false; // TODO:
        if (concordance) {
            filters.add("HMF_ClippingConcordance");
        }

        return filters.isEmpty() ? "PASS" : String.join(";", filters);
    }

}
