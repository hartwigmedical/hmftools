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
import com.hartwig.hmftools.breakpointinspector.datamodel.EnrichedVariantContext;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.vcf.VCFHeader;

final class Filter {

    private static final int MIN_ANCHOR_LENGTH = 30;

    private Filter() {
    }

    static void updateVCFHeader(@NotNull final VCFHeader header) {
        Arrays.stream(FilterType.values()).forEach(filterType -> header.addMetaDataLine(filterType.toHeaderLine()));
    }

    @NotNull
    static Collection<String> errorFilter() {
        return Collections.singletonList(FilterType.BREAKPOINT_ERROR.toString());
    }

    @NotNull
    static Collection<String> filters(final EnrichedVariantContext variant, final SampleStats tumorStats, final SampleStats refStats,
            final Pair<Location, Location> breakpoints, final float contamination) {
        final List<FilterType> filters = Lists.newArrayList();

        if (Stream.of(tumorStats.BP1_Stats, tumorStats.BP2_Stats)
                .mapToInt(s -> s.PR_Only_Normal + s.PR_SR_Normal + s.PR_Only_Support + s.PR_SR_Support)
                .anyMatch(i -> i < 10)) {
            filters.add(FilterType.MIN_DEPTH);
        }

        final int tumorSR = Stream.of(tumorStats.BP1_Stats, tumorStats.BP2_Stats).mapToInt(Filter::supportSR).sum();

        if (variant.isShortVariant()) {
            final boolean bothSidesHaveSR = Stream.of(tumorStats.BP1_Stats, tumorStats.BP2_Stats).allMatch(s -> supportSR(s) > 0);
            final boolean anchorLengthOkay = tumorStats.SR_Evidence.stream()
                    .anyMatch(p -> Stream.of(p.getLeft(), p.getRight())
                            .anyMatch(r -> r.getAlignmentEnd() - r.getAlignmentStart() >= MIN_ANCHOR_LENGTH));

            if (!bothSidesHaveSR) {
                filters.add(FilterType.SR_SUPPORT_ZERO);
            } else if (!anchorLengthOkay) {
                filters.add(FilterType.MIN_ANCHOR_LENGTH);
            }

            // NERA: must not have SR support in normal
            final int refSR = Stream.of(refStats.BP1_Stats, refStats.BP2_Stats).mapToInt(Filter::supportSR).sum();
            final int allowableNormalSupport = (int) (contamination * tumorSR);
            if (refSR > allowableNormalSupport) {
                filters.add(FilterType.SR_NORMAL_SUPPORT);
            }
        } else if (variant.isInsert()) {
            // NERA: we only need to check BP1 as BP1 PR+PRSR == BP2 PR+PRSR
            final int allowableNormalSupport = (int) (contamination * supportPR(tumorStats.BP1_Stats));
            if (supportPR(refStats.BP1_Stats) > allowableNormalSupport) {
                filters.add(FilterType.PR_NORMAL_SUPPORT);
            }

            final boolean anchorLengthOkay = tumorStats.PR_Evidence.stream()
                    .anyMatch(p -> Stream.of(p.getLeft(), p.getRight())
                            .allMatch(r -> r.getAlignmentEnd() - r.getAlignmentStart() >= MIN_ANCHOR_LENGTH));

            // NERA: only applicable for longer variants
            final int tumorPR = Stream.of(tumorStats.BP1_Stats, tumorStats.BP2_Stats).mapToInt(Filter::supportPR).sum();
            if (tumorPR == 0) {
                filters.add(FilterType.PR_SUPPORT_ZERO);
            } else if (!anchorLengthOkay) {
                filters.add(FilterType.MIN_ANCHOR_LENGTH);
            }
        }

        // NERA: We must adjust from Manta breakpoint convention to our clipping position convention
        final List<Location> adjustedBP =
                Arrays.asList(breakpoints.getLeft().add(variant.orientationBP1()), breakpoints.getRight().add(variant.orientationBP2()));

        final Set<String> concordantReads = Sets.newHashSet();
        for (final Location bp : adjustedBP) {
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
                            concordantReads.addAll(r.SupportingReads);
                        }
                    } else {
                        if (tumorSeq.equals(r.LongestClipSequence.substring(0, 5))) {
                            concordantReads.addAll(r.SupportingReads);
                        }
                    }
                }
            }
        }

        if (concordantReads.size() > (int) (contamination * tumorSR)) {
            filters.add(FilterType.CLIPPING_CONCORDANCE);
        }

        final Set<String> merged = Sets.newHashSet(variant.filters());
        merged.addAll(filters.stream().map(FilterType::toString).collect(Collectors.toList()));
        return merged;
    }

    private static int supportPR(@NotNull final BreakpointStats stats) {
        return stats.PR_SR_Support + stats.PR_Only_Support;
    }

    private static int supportSR(@NotNull final BreakpointStats stats) {
        return stats.PR_SR_Support + stats.SR_Only_Support;
    }
}
