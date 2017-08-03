package com.hartwig.hmftools.breakpointinspector;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;

import org.apache.commons.lang3.tuple.Pair;

import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

class Filter {

    private enum Filters {
        BreakpointError("HMF_BreakpointError", "BPI failed to determine breakpoints"),
        MinDepth("HMF_MinDepth", "The depth across one of the breakpoints is <10"),
        MinAnchorLength("HMF_MinAnchorLength", "There isn't at least one PR with >=30 bases matched in alignment"),
        SRSupportZero("HMF_SRSupportZero", "Short delete (<2000) must have SR support"),
        SRNormalSupport("HMF_SRNormalSupport", "Short delete (<2000) has SR support in normal"),
        PRNormalSupport("HMF_PRNormalSupport", "PR support in the normal"),
        PRSupportZero("HMF_PRSupportZero", "No PR support in tumor"),
        ClippingConcordance("HMF_ClippingConcordance", "At least 5 base clipped bases concordance between tumor and normal");

        private final String Name;
        private final String Description;

        Filters(final String name, final String description) {
            Name = name;
            Description = description;
        }

        VCFFilterHeaderLine toHeaderLine() {
            return new VCFFilterHeaderLine(Name, Description);
        }

        @Override
        public String toString() {
            return Name;
        }
    }

    static void updateHeader(final VCFHeader header) {
        Arrays.stream(Filters.values()).forEach(f -> header.addMetaDataLine(f.toHeaderLine()));
    }

    static Collection<String> getErrorFilter() {
        return Collections.singletonList(Filters.BreakpointError.toString());
    }

    static Collection<String> getFilters(final HMFVariantContext ctx, final SampleStats tumorStats, final SampleStats refStats,
            final Pair<Location, Location> breakpoints) {

        final List<Filters> filters = Lists.newArrayList();

        if (Stream.of(tumorStats.BP1_Stats, tumorStats.BP2_Stats)
                .mapToInt(s -> s.PR_Only_Normal + s.PR_SR_Normal + s.PR_Only_Support + s.PR_SR_Support)
                .anyMatch(i -> i < 10)) {
            filters.add(Filters.MinDepth);
        }

        final boolean anchorLengthOkay = tumorStats.PR_Evidence.stream()
                .anyMatch(p -> Stream.of(p.getLeft(), p.getRight()).anyMatch(r -> r.getAlignmentEnd() - r.getAlignmentStart() >= 30));
        if (!anchorLengthOkay) {
            filters.add(Filters.MinAnchorLength);
        }

        if (ctx.isShortDelete()) {
            // short delete logic, must have SR support
            final int tumor_SR =
                    Stream.of(tumorStats.BP1_Stats, tumorStats.BP2_Stats).mapToInt(s -> s.PR_SR_Support + s.SR_Only_Support).sum();
            if (tumor_SR == 0) {
                filters.add(Filters.SRSupportZero);
            }

            // short delete logic, must not have SR support in normal
            final int ref_SR = Stream.of(refStats.BP1_Stats, refStats.BP2_Stats).mapToInt(s -> s.PR_SR_Support + s.SR_Only_Support).sum();
            if (ref_SR > 0) {
                filters.add(Filters.SRNormalSupport);
            }
        } else {
            // we only need to check BP1 as BP1 PR+PRSR == BP2 PR+PRSR
            if (refStats.BP1_Stats.PR_Only_Support + refStats.BP1_Stats.PR_SR_Support > 0) {
                filters.add(Filters.PRNormalSupport);
            }
        }

        if (Stream.of(tumorStats.BP1_Stats, tumorStats.BP2_Stats).mapToInt(s -> s.PR_Only_Support + s.PR_SR_Support).sum() == 0) {
            filters.add(Filters.PRSupportZero);
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
            filters.add(Filters.ClippingConcordance);
        }

        final List<String> merged = Lists.newArrayList(ctx.Filter);
        merged.addAll(filters.stream().map(Filters::toString).collect(Collectors.toList()));
        return merged;
    }

}
