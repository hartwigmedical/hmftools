package com.hartwig.hmftools.breakpointinspector;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Stream;

class Filter {

    private static boolean concordantStrings(final String a, final String b) {
        if (a == null || b == null)
            return false;
        if (a.length() < b.length())
            return !a.isEmpty() && b.startsWith(a) || b.endsWith(a);
        else
            return !b.isEmpty() && a.startsWith(b) || a.endsWith(b);
    }

    static String getFilterString(final Util.HMFVariantContext ctx, final Stats.SampleStats tumorStats,
            final Stats.SampleStats refStats) {

        final List<String> filters = new ArrayList<>(ctx.Filter);

        if (ctx.Type == Util.HMFVariantType.DEL && ctx.MantaBP1.ReferenceIndex == ctx.MantaBP2.ReferenceIndex
                && (ctx.MantaBP2.Position - ctx.MantaBP1.Position) < 2000) {
            // short delete logic, must have SR support
            final int SR = Stream.of(tumorStats.BP1_Stats, tumorStats.BP2_Stats).mapToInt(
                    s -> s.PR_SR_Support + s.SR_Only_Support).sum();
            if (SR == 0) {
                filters.add("HMF_SRSupportZero");
            }
        } else {
            // we only need to check BP1 as BP1 PR+PRSR == BP2 PR+PRSR
            if (refStats.BP1_Stats.PR_Only_Support + refStats.BP1_Stats.PR_SR_Support > 0) {
                filters.add("HMF_PRNormalSupport");
            }
        }

        boolean concordance = false;
        for (final Util.Location bp : Arrays.asList(tumorStats.BP1, tumorStats.BP2)) {
            final Stats.Clip tumor_clip = tumorStats.Clipping_Stats.LocationMap.get(bp);
            final Stats.Clip ref_clip = refStats.Clipping_Stats.LocationMap.get(bp);
            concordance |= tumor_clip != null && ref_clip != null && concordantStrings(tumor_clip.LongestClipSequence,
                    ref_clip.LongestClipSequence);
        }
        if (concordance)
            filters.add("HMF_ClippingConcordance");

        return filters.isEmpty() ? "PASS" : String.join(";", filters);
    }

}
