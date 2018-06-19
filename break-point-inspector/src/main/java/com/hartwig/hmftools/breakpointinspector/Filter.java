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
import com.hartwig.hmftools.breakpointinspector.datamodel.BreakpointStats;
import com.hartwig.hmftools.breakpointinspector.datamodel.EnrichedVariantContext;
import com.hartwig.hmftools.breakpointinspector.datamodel.SampleStats;

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
            final Pair<Location, Location> breakpoints, final double contaminationFraction) {
        final List<FilterType> filters = Lists.newArrayList();

        if (Stream.of(tumorStats.bp1Stats, tumorStats.bp2Stats)
                .mapToInt(s -> s.prOnlyNormal + s.prSrNormal + s.prOnlySupport + s.prSrSupport)
                .anyMatch(i -> i < 10)) {
            filters.add(FilterType.MIN_DEPTH);
        }

        final int tumorSR = Stream.of(tumorStats.bp1Stats, tumorStats.bp2Stats).mapToInt(Filter::supportSR).sum();

        if (variant.isShortVariant()) {
            final boolean bothSidesHaveSR = Stream.of(tumorStats.bp1Stats, tumorStats.bp2Stats).allMatch(s -> supportSR(s) > 0);
            final boolean anchorLengthOkay = tumorStats.srEvidence.stream()
                    .anyMatch(p -> Stream.of(p.getLeft(), p.getRight())
                            .anyMatch(r -> r.getAlignmentEnd() - r.getAlignmentStart() >= MIN_ANCHOR_LENGTH));

            if (!bothSidesHaveSR) {
                filters.add(FilterType.SR_SUPPORT_ZERO);
            } else if (!anchorLengthOkay) {
                filters.add(FilterType.MIN_ANCHOR_LENGTH);
            }

            // NERA: must not have SR support in normal
            final int refSR = Stream.of(refStats.bp1Stats, refStats.bp2Stats).mapToInt(Filter::supportSR).sum();
            final int allowableNormalSupport = (int) (contaminationFraction * tumorSR);
            if (refSR > allowableNormalSupport) {
                filters.add(FilterType.SR_NORMAL_SUPPORT);
            }
        } else if (!variant.isInsert()) {
            // NERA: we only need to check BP1 as BP1 PR+PRSR == BP2 PR+PRSR
            final int allowableNormalSupport = (int) (contaminationFraction * supportPR(tumorStats.bp1Stats));
            if (supportPR(refStats.bp1Stats) > allowableNormalSupport) {
                filters.add(FilterType.PR_NORMAL_SUPPORT);
            }

            final boolean anchorLengthOkay = tumorStats.prEvidence.stream()
                    .anyMatch(p -> Stream.of(p.getLeft(), p.getRight())
                            .allMatch(r -> r.getAlignmentEnd() - r.getAlignmentStart() >= MIN_ANCHOR_LENGTH));

            // NERA: only applicable for longer variants
            final int tumorPR = Stream.of(tumorStats.bp1Stats, tumorStats.bp2Stats).mapToInt(Filter::supportPR).sum();
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
            for (final ClipStats t : tumorStats.sampleClipping.sequencesAt(bp)) {
                if (t.longestClipSequence.length() < 5) {
                    continue;
                }

                final String tumorSeq = t.left
                        ? t.longestClipSequence.substring(t.longestClipSequence.length() - 5)
                        : t.longestClipSequence.substring(0, 5);

                for (final ClipStats r : refStats.sampleClipping.sequencesAt(bp)) {
                    if (t.left != r.left) {
                        continue;
                    } else if (r.longestClipSequence.length() < 5) {
                        continue;
                    }

                    if (t.left) {
                        if (tumorSeq.equals(r.longestClipSequence.substring(r.longestClipSequence.length() - 5))) {
                            concordantReads.addAll(r.supportingReads);
                        }
                    } else {
                        if (tumorSeq.equals(r.longestClipSequence.substring(0, 5))) {
                            concordantReads.addAll(r.supportingReads);
                        }
                    }
                }
            }
        }

        if (concordantReads.size() > (int) (contaminationFraction * tumorSR)) {
            filters.add(FilterType.CLIPPING_CONCORDANCE);
        }

        final Set<String> merged = Sets.newHashSet(variant.filters());
        merged.addAll(filters.stream().map(FilterType::toString).collect(Collectors.toList()));
        return merged;
    }

    private static int supportPR(@NotNull final BreakpointStats stats) {
        return stats.prSrSupport + stats.prOnlySupport;
    }

    private static int supportSR(@NotNull final BreakpointStats stats) {
        return stats.prSrSupport + stats.srOnlySupport;
    }
}
