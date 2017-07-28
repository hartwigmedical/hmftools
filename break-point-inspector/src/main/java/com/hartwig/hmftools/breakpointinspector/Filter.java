package com.hartwig.hmftools.breakpointinspector;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Stream;

import com.google.common.collect.Lists;

import org.apache.commons.lang3.tuple.Pair;

class Filter {

    static String getFilterString(final HMFVariantContext ctx, final SampleStats tumorStats, final SampleStats refStats,
            final Pair<Location, Location> breakpoints) {

        final List<String> filters = Lists.newArrayList(ctx.Filter);

        if (Stream.of(tumorStats.BP1_Stats, tumorStats.BP2_Stats)
                .mapToInt(s -> s.PR_Only_Normal + s.PR_SR_Normal + s.PR_Only_Support + s.PR_SR_Support)
                .anyMatch(i -> i < 10)) {
            filters.add("HMF_MinDepth");
        }

        final boolean anchorLengthOkay = tumorStats.PR_Evidence.stream()
                .anyMatch(p -> Stream.of(p.getLeft(), p.getRight()).anyMatch(r -> r.getAlignmentEnd() - r.getAlignmentStart() >= 30));
        if (!anchorLengthOkay) {
            filters.add("HMF_MinAnchorLength");
        }

        if (ctx.Type == HMFVariantType.DEL && ctx.MantaBP1.ReferenceIndex == ctx.MantaBP2.ReferenceIndex
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

        if (Stream.of(tumorStats.BP1_Stats, tumorStats.BP2_Stats).mapToInt(s -> s.PR_Only_Support + s.PR_SR_Support).sum() == 0) {
            filters.add("HMF_PRSupportZero");
        }

        final List<Location> adjusted_bp = Arrays.asList(breakpoints.getLeft().add(ctx.OrientationBP1 > 0 ? 1 : 0),
                breakpoints.getRight().add(ctx.OrientationBP2 > 0 ? 1 : 0));

        boolean concordance = false;
        for (final Location bp : adjusted_bp) {
            for (final ClipStats t : tumorStats.Sample_Clipping.getSequencesAt(bp)) {

                if (t.LongestClipSequence.length() < 5) {
                    continue;
                }

                final String tumorSeq = t.Left
                        ? t.LongestClipSequence.substring(t.LongestClipSequence.length() - 5)
                        : t.LongestClipSequence.substring(0, 5);

                for (final ClipStats r : refStats.Sample_Clipping.getSequencesAt(bp)) {

                    if (t.Left != r.Left) {
                        continue;
                    } else if (r.LongestClipSequence.length() < 5) {
                        continue;
                    }

                    if (t.Left) {
                        concordance |= tumorSeq.equals(r.LongestClipSequence.substring(r.LongestClipSequence.length() - 5));
                    } else {
                        concordance |= tumorSeq.equals(r.LongestClipSequence.substring(0, 5));
                    }
                }
            }
        }

        if (concordance) {
            filters.add("HMF_ClippingConcordance");
        }

        return filters.isEmpty() ? "PASS" : String.join(";", filters);
    }

}
