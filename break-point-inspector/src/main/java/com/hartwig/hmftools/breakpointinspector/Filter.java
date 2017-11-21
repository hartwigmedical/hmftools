package com.hartwig.hmftools.breakpointinspector;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.breakpointinspector.clipping.ClipStats;

import org.apache.commons.lang3.tuple.Pair;

import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

class Filter {

    enum Filters {
        BreakpointError("BPI_BreakpointError", "BPI failed to determine breakpoints"),
        MinDepth("BPI_MinDepth", "The depth across one of the breakpoints is <10"),
        MinAnchorLength("BPI_MinAnchorLength", "There isn't at least one PR with >=30 bases matched in both alignments"),
        SRSupportZero("BPI_SRSupportZero", "Short delete or dupe (<1000) must have SR support"),
        SRNormalSupport("BPI_SRNormalSupport", "Short delete or dupe (<1000) has SR support in normal"),
        PRNormalSupport("BPI_PRNormalSupport", "PR support in the normal"),
        PRSupportZero("BPI_PRSupportZero", "No PR support in tumor"),
        ClippingConcordance("BPI_ClippingConcordance", "At least 5 base clipped bases concordance between tumor and normal");

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

    static void UpdateVCFHeader(final VCFHeader header) {
        Arrays.stream(Filters.values()).forEach(f -> header.addMetaDataLine(f.toHeaderLine()));
    }

    static Collection<String> getErrorFilter() {
        return Collections.singletonList(Filters.BreakpointError.toString());
    }

    private static int supportPR(final BreakpointStats s) {
        return s.PR_SR_Support + s.PR_Only_Support;
    }

    private static int supportSR(final BreakpointStats s) {
        return s.PR_SR_Support + s.SR_Only_Support;
    }

    static Collection<String> getFilters(final HMFVariantContext ctx, final SampleStats tumorStats, final SampleStats refStats,
            final Pair<Location, Location> breakpoints, final float contamination) {

        final int MIN_ANCHOR_LENGTH = 30;

        final List<Filters> filters = Lists.newArrayList();

        if (Stream.of(tumorStats.BP1_Stats, tumorStats.BP2_Stats)
                .mapToInt(s -> s.PR_Only_Normal + s.PR_SR_Normal + s.PR_Only_Support + s.PR_SR_Support)
                .anyMatch(i -> i < 10)) {
            filters.add(Filters.MinDepth);
        }

        final int tumor_SR = Stream.of(tumorStats.BP1_Stats, tumorStats.BP2_Stats).mapToInt(Filter::supportSR).sum();

        if (ctx.isInsert()) {

            // no PR/SR checks

        } else if (ctx.isShortDelete() || ctx.isShortDuplicate()) {
            // short variant logic

            final boolean bothSidesHaveSR = Stream.of(tumorStats.BP1_Stats, tumorStats.BP2_Stats).allMatch(s -> supportSR(s) > 0);
            final boolean anchorLengthOkay = tumorStats.SR_Evidence.stream()
                    .anyMatch(p -> Stream.of(p.getLeft(), p.getRight())
                            .anyMatch(r -> r.getAlignmentEnd() - r.getAlignmentStart() >= MIN_ANCHOR_LENGTH));

            if (!bothSidesHaveSR) {
                filters.add(Filters.SRSupportZero);
            } else if (!anchorLengthOkay) {
                filters.add(Filters.MinAnchorLength);
            }

            // must not have SR support in normal
            final int ref_SR = Stream.of(refStats.BP1_Stats, refStats.BP2_Stats).mapToInt(Filter::supportSR).sum();
            final int allowableNormalSupport = (int) (contamination * tumor_SR);
            if (ref_SR > allowableNormalSupport) {
                filters.add(Filters.SRNormalSupport);
            }
        } else {

            // we only need to check BP1 as BP1 PR+PRSR == BP2 PR+PRSR
            final int allowableNormalSupport = (int) (contamination * supportPR(tumorStats.BP1_Stats));
            if (supportPR(refStats.BP1_Stats) > allowableNormalSupport) {
                filters.add(Filters.PRNormalSupport);
            }

            final boolean anchorLengthOkay = tumorStats.PR_Evidence.stream()
                    .anyMatch(p -> Stream.of(p.getLeft(), p.getRight())
                            .allMatch(r -> r.getAlignmentEnd() - r.getAlignmentStart() >= MIN_ANCHOR_LENGTH));

            // only applicable for longer variants
            final int tumor_PR = Stream.of(tumorStats.BP1_Stats, tumorStats.BP2_Stats).mapToInt(Filter::supportPR).sum();
            if (tumor_PR == 0) {
                filters.add(Filters.PRSupportZero);
            } else if (!anchorLengthOkay) {
                filters.add(Filters.MinAnchorLength);
            }
        }

        // we must adjust from Manta breakpoint convention to our clipping position convention
        final List<Location> adjusted_bp =
                Arrays.asList(breakpoints.getLeft().add(ctx.OrientationBP1), breakpoints.getRight().add(ctx.OrientationBP2));

        final Set<String> concordant_reads = Sets.newHashSet();
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
                        if (tumorSeq.equals(r.LongestClipSequence.substring(r.LongestClipSequence.length() - 5))) {
                            concordant_reads.addAll(r.SupportingReads);
                        }
                    } else {
                        if (tumorSeq.equals(r.LongestClipSequence.substring(0, 5))) {
                            concordant_reads.addAll(r.SupportingReads);
                        }
                    }
                }

            }

        }

        if (concordant_reads.size() > (int) (contamination * tumor_SR)) {
            filters.add(Filters.ClippingConcordance);
        }

        final Set<String> merged = Sets.newHashSet(ctx.Filter);
        merged.addAll(filters.stream().map(Filters::toString).collect(Collectors.toList()));
        return merged;
    }

}
