package com.hartwig.hmftools.breakpointinspector;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.breakpointinspector.clipping.Clipping;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

class Analysis {

    private static final Logger LOGGER = LogManager.getLogger(Analysis.class);

    private final SamReader refReader;
    private final SamReader tumorReader;

    @Nullable
    private final SAMFileWriter refWriter;
    @Nullable
    private final SAMFileWriter tumorWriter;

    private final int range;
    private final int[] extraUncertainties;

    private final Set<Integer> tumorWrittenReads = new HashSet<>();
    private final Set<Integer> refWrittenReads = new HashSet<>();

    Analysis(final SamReader refReader, @Nullable final SAMFileWriter refWriter, final SamReader tumorReader,
            @Nullable final SAMFileWriter tumorWriter, final int range, final int[] extraUncertainties) {
        this.refReader = refReader;
        this.refWriter = refWriter;
        this.tumorReader = tumorReader;
        this.tumorWriter = tumorWriter;
        this.range = range;
        this.extraUncertainties = extraUncertainties;
    }

    private static class AlignmentList extends ArrayList<SAMRecord> {
    }

    private static class AlignmentMap extends HashMap<String, AlignmentList> {
        AlignmentMap(int initialSize) {
            super(initialSize);
        }
    }

    private static class PairedReads extends ArrayList<Pair<SAMRecord, SAMRecord>> {
        PairedReads(int initialSize) {
            super(initialSize);
        }
    }

    private static int orientation(final SAMRecord record) {
        return record.getReadNegativeStrandFlag() ? -1 : 1;
    }

    private static Pair<Integer, Integer> orientation(final Pair<SAMRecord, SAMRecord> pair) {
        return Pair.of(orientation(pair.getLeft()), orientation(pair.getRight()));
    }

    private static <L> Stream<L> stream(final Pair<L, L> pair) {
        return Stream.of(pair.getLeft(), pair.getRight());
    }

    private static boolean isMate(final SAMRecord read, final SAMRecord mate) {
        return read.getReadName().equals(mate.getReadName()) && read.getMateReferenceIndex().equals(mate.getReferenceIndex())
                && Math.abs(read.getMateAlignmentStart() - mate.getAlignmentStart()) <= 1;
    }

    private static boolean span(final Pair<SAMRecord, SAMRecord> pair, final Location breakpoint) {
        return Location.fromSAMRecord(pair.getLeft(), true).compareTo(breakpoint) <= 0
                && Location.fromSAMRecord(pair.getRight(), false).compareTo(breakpoint) >= 0;
    }

    private static boolean overlap(final SAMRecord read, final Location breakpoint) {
        return read.getReferenceIndex() == breakpoint.ReferenceIndex && read.getAlignmentStart() <= breakpoint.Position
                && breakpoint.Position <= read.getAlignmentEnd();
    }

    private static boolean clippedOnCorrectSide(final SAMRecord record, final int orientation) {
        return (orientation > 0 ? record.getCigar().isRightClipped() : record.getCigar().isLeftClipped());
    }

    private boolean withinRange(final Location a, final Location b, final Range range, final int extraUncertainty) {
        return a.ReferenceIndex == b.ReferenceIndex && (a.Position >= b.Position + range.Start - extraUncertainty) && (a.Position
                <= b.Position + range.End + extraUncertainty);
    }

    private SampleStats collectEvidence(final HMFVariantContext ctx, final PairedReads pairs, final Pair<Location, Location> breakpoints) {

        final SampleStats result = new SampleStats();
        final Pair<Integer, Integer> ctxOrientation = Pair.of(ctx.OrientationBP1, ctx.OrientationBP2);

        for (final Pair<SAMRecord, SAMRecord> pair : pairs) {

            final boolean proper = stream(pair).anyMatch(SAMRecord::getProperPairFlag);
            final boolean secondary = stream(pair).anyMatch(SAMRecord::isSecondaryOrSupplementary);
            final boolean correctOrientation = orientation(pair).equals(ctxOrientation);
            final boolean correctChromosome =
                    Location.fromSAMRecord(pair.getLeft()).sameChromosomeAs(ctx.MantaBP1) && Location.fromSAMRecord(pair.getRight())
                            .sameChromosomeAs(ctx.MantaBP2);

            final int breakProximity = 200;
            final boolean leftCorrectPosition = ctx.OrientationBP1 > 0
                    ? breakpoints.getLeft().Position - pair.getLeft().getAlignmentEnd() < breakProximity
                    : pair.getLeft().getAlignmentStart() - breakpoints.getLeft().Position < breakProximity;
            final boolean rightCorrectPosition = ctx.OrientationBP2 > 0
                    ? breakpoints.getRight().Position - pair.getRight().getAlignmentEnd() < breakProximity
                    : pair.getRight().getAlignmentStart() - breakpoints.getRight().Position < breakProximity;

            boolean support = correctOrientation && correctChromosome && leftCorrectPosition && rightCorrectPosition;
            if (support) {

                final int left_outer = Location.fromSAMRecord(pair.getLeft(), ctx.OrientationBP1 > 0).compareTo(breakpoints.getLeft());
                final int left_inner = Location.fromSAMRecord(pair.getLeft(), ctx.OrientationBP1 < 0).compareTo(breakpoints.getLeft());
                final int right_outer = Location.fromSAMRecord(pair.getRight(), ctx.OrientationBP2 > 0).compareTo(breakpoints.getRight());
                final int right_inner = Location.fromSAMRecord(pair.getRight(), ctx.OrientationBP2 < 0).compareTo(breakpoints.getRight());

                if (ctx.OrientationBP1 > 0) {
                    support &= left_outer < 0;
                    support &= left_inner <= 0;
                } else {
                    support &= left_outer > 0;
                    support &= left_inner >= 0;
                }
                if (ctx.OrientationBP2 > 0) {
                    support &= right_outer < 0;
                    support &= right_inner <= 0;
                } else {
                    support &= right_outer > 0;
                    support &= right_inner >= 0;
                }

            }

            if (support) {

                final boolean clip_bp1 = Location.fromSAMRecord(pair.getLeft(), ctx.OrientationBP1 < 0)
                        .add(ctx.OrientationBP1 < 0 ? -1 : 0)
                        .compareTo(breakpoints.getLeft()) == 0 && clippedOnCorrectSide(pair.getLeft(), ctx.OrientationBP1);

                final boolean clip_bp2 = Location.fromSAMRecord(pair.getRight(), ctx.OrientationBP2 < 0)
                        .add(ctx.OrientationBP2 < 0 ? -1 : 0)
                        .compareTo(breakpoints.getRight()) == 0 && clippedOnCorrectSide(pair.getRight(), ctx.OrientationBP2);

                if (clip_bp1) {
                    result.BP1_Stats.PR_SR_Support++;
                } else {
                    result.BP1_Stats.PR_Only_Support++;
                }

                if (clip_bp2) {
                    result.BP2_Stats.PR_SR_Support++;
                } else {
                    result.BP2_Stats.PR_Only_Support++;
                }

                result.PR_Evidence.add(pair);

            } else if (proper || secondary) {

                final boolean clip_bp1 =
                        Location.fromSAMRecord(ctx.OrientationBP1 > 0 ? pair.getRight() : pair.getLeft(), ctx.OrientationBP1 < 0)
                                .add(ctx.OrientationBP1 > 0 ? 0 : -1)
                                .compareTo(breakpoints.getLeft()) == 0 && clippedOnCorrectSide(
                                ctx.OrientationBP1 > 0 ? pair.getRight() : pair.getLeft(), ctx.OrientationBP1);
                final boolean clip_bp2 =
                        Location.fromSAMRecord(ctx.OrientationBP2 > 0 ? pair.getRight() : pair.getLeft(), ctx.OrientationBP2 < 0)
                                .add(ctx.OrientationBP2 > 0 ? 0 : -1)
                                .compareTo(breakpoints.getRight()) == 0 && clippedOnCorrectSide(
                                ctx.OrientationBP2 > 0 ? pair.getRight() : pair.getLeft(), ctx.OrientationBP2);

                final boolean span_bp1 = span(pair, breakpoints.getLeft());
                final boolean span_bp2 = span(pair, breakpoints.getRight());
                final boolean overlap_bp1 = stream(pair).anyMatch(r -> overlap(r, breakpoints.getLeft()));
                final boolean overlap_bp2 = stream(pair).anyMatch(r -> overlap(r, breakpoints.getRight()));

                if (clip_bp1) {
                    result.BP1_Stats.SR_Only_Support++;
                    continue;
                }
                if (clip_bp2) {
                    result.BP2_Stats.SR_Only_Support++;
                    continue;
                }

                if (span_bp1) {
                    if (overlap_bp1) {
                        result.BP1_Stats.PR_SR_Normal++;
                    } else {
                        result.BP1_Stats.PR_Only_Normal++;
                    }
                }
                if (span_bp2) {
                    if (overlap_bp2) {
                        result.BP2_Stats.PR_SR_Normal++;
                    } else {
                        result.BP2_Stats.PR_Only_Normal++;
                    }
                }

            }
        }

        return result;
    }

    enum BreakpointError {
        NONE,
        ALGO_ERROR,
        NO_EVIDENCE,
        UNINITIALIZED
    }

    private static class BreakpointResult {
        private BreakpointResult(final Pair<Location, Location> breakpoints) {
            Breakpoints = breakpoints;
        }

        private BreakpointResult(final BreakpointError error) {
            Breakpoints = Pair.of(null, null);
            Error = error;
        }

        static BreakpointResult from(final Pair<Location, Location> breakpoints) {
            return new BreakpointResult(breakpoints);
        }

        static BreakpointResult from(final BreakpointError error) {
            return new BreakpointResult(error);
        }

        Pair<Location, Location> Breakpoints;
        BreakpointError Error = BreakpointError.NONE;
    }

    private BreakpointResult determineBreakpoints(final HMFVariantContext ctx, final PairedReads pairs, final int extraUncertainty) {

        final Pair<Integer, Integer> ctxOrientation = Pair.of(ctx.OrientationBP1, ctx.OrientationBP2);

        // extract all interesting pairs

        final PairedReads interesting = new PairedReads(pairs.size() / 2);
        final PairedReads clipped_proper = new PairedReads(pairs.size() / 2);
        final PairedReads secondary_pairs = new PairedReads(pairs.size() / 2);

        for (final Pair<SAMRecord, SAMRecord> pair : pairs) {

            final boolean correctOrientation = orientation(pair).equals(ctxOrientation);
            final boolean correctChromosome =
                    Location.fromSAMRecord(pair.getLeft()).sameChromosomeAs(ctx.MantaBP1) && Location.fromSAMRecord(pair.getRight())
                            .sameChromosomeAs(ctx.MantaBP2);
            final boolean hasExpectedClipping =
                    clippedOnCorrectSide(pair.getLeft(), ctx.OrientationBP1) || clippedOnCorrectSide(pair.getRight(), ctx.OrientationBP2);

            final boolean proper = stream(pair).anyMatch(SAMRecord::getProperPairFlag);
            final boolean potentialSROnly = Stream.of(ctx.OrientationBP1, ctx.OrientationBP2)
                    .anyMatch(orientation -> clippedOnCorrectSide(orientation > 0 ? pair.getRight() : pair.getLeft(), orientation));

            final boolean secondary = stream(pair).anyMatch(SAMRecord::isSecondaryOrSupplementary);

            LOGGER.debug(
                    "{} {}->{} {} {}->{} correctOrientation({}) correctChromosome({}) hasExpectedClipping({}) proper({}) potentialSROnly({}) secondary({})",
                    pair.getLeft().getReadName(), pair.getLeft().getAlignmentStart(), pair.getLeft().getMateAlignmentStart(),
                    pair.getRight().getReadName(), pair.getRight().getAlignmentStart(), pair.getRight().getMateAlignmentStart(),
                    correctOrientation, correctChromosome, hasExpectedClipping, proper, potentialSROnly, secondary);

            // TODO: check insert size?

            if (secondary && potentialSROnly) {
                secondary_pairs.add(pair);
            } else if ((!proper || hasExpectedClipping) && correctChromosome && correctOrientation) {
                interesting.add(pair);
            } else if (proper && potentialSROnly) {
                clipped_proper.add(pair);
            }

        }

        // load clipping info

        Clipping bp1_clipping = new Clipping();
        Clipping bp2_clipping = new Clipping();

        for (final Pair<SAMRecord, SAMRecord> pair : interesting) {
            if (ctx.OrientationBP1 > 0) {
                bp1_clipping.add(Clipping.getRightClip(pair.getLeft()));
            } else {
                bp1_clipping.add(Clipping.getLeftClip(pair.getLeft()));
            }
            if (ctx.OrientationBP2 > 0) {
                bp2_clipping.add(Clipping.getRightClip(pair.getRight()));
            } else {
                bp2_clipping.add(Clipping.getLeftClip(pair.getRight()));
            }
        }

        // include more clipping information

        for (final Pair<SAMRecord, SAMRecord> pair : clipped_proper) {
            if (stream(pair).allMatch(r -> Location.fromSAMRecord(r).sameChromosomeAs(ctx.MantaBP1))) {
                if (ctx.OrientationBP1 > 0) {
                    bp1_clipping.add(Clipping.getRightClip(pair.getRight()));
                } else {
                    bp1_clipping.add(Clipping.getLeftClip(pair.getLeft()));
                }
            }
            if (stream(pair).allMatch(r -> Location.fromSAMRecord(r).sameChromosomeAs(ctx.MantaBP2))) {
                if (ctx.OrientationBP2 > 0) {
                    bp2_clipping.add(Clipping.getRightClip(pair.getRight()));
                } else {
                    bp2_clipping.add(Clipping.getLeftClip(pair.getLeft()));
                }
            }
        }

        // include secondary clipping information

        for (final Pair<SAMRecord, SAMRecord> pair : secondary_pairs) {
            if (stream(pair).allMatch(r -> Location.fromSAMRecord(r).sameChromosomeAs(ctx.MantaBP1))) {
                if (ctx.OrientationBP1 > 0) {
                    bp1_clipping.add(Clipping.getRightClip(pair.getRight()));
                } else {
                    bp1_clipping.add(Clipping.getLeftClip(pair.getLeft()));
                }
            }
            if (stream(pair).allMatch(r -> Location.fromSAMRecord(r).sameChromosomeAs(ctx.MantaBP2))) {
                if (ctx.OrientationBP2 > 0) {
                    bp2_clipping.add(Clipping.getRightClip(pair.getRight()));
                } else {
                    bp2_clipping.add(Clipping.getLeftClip(pair.getLeft()));
                }
            }
        }

        // determinate candidates based on clipping info

        final List<Location> bp1_candidates = bp1_clipping.getSequences()
                .stream()
                //.filter(c -> c.LongestClipSequence.length() >= 5)
                .map(c -> c.Alignment)
                .filter(c -> withinRange(c, ctx.MantaBP1, ctx.Uncertainty1, extraUncertainty))
                .collect(Collectors.toList());

        if (bp1_candidates.isEmpty()) {
            final Location candidate = interesting.stream()
                    .map(Pair::getLeft)
                    .map(r -> Location.fromSAMRecord(r, ctx.OrientationBP1 < 0).add(ctx.OrientationBP1 > 0 ? 0 : -1))
                    .filter(l -> withinRange(l, ctx.MantaBP1, ctx.Uncertainty1, extraUncertainty))
                    .max((a, b) -> ctx.OrientationBP1 > 0 ? a.compareTo(b) : b.compareTo(a))
                    .orElse(null);

            if (candidate == null) {
                return BreakpointResult.from(BreakpointError.NO_EVIDENCE);
            }

            bp1_candidates.add(candidate);
        }

        final List<Location> bp2_candidates = bp2_clipping.getSequences()
                .stream()
                //.filter(c -> c.LongestClipSequence.length() >= 5)
                .map(c -> c.Alignment)
                .filter(c -> withinRange(c, ctx.MantaBP2, ctx.Uncertainty2, extraUncertainty))
                .collect(Collectors.toList());

        if (bp2_candidates.isEmpty()) {
            final Location candidate = interesting.stream()
                    .map(Pair::getRight)
                    .map(r -> Location.fromSAMRecord(r, ctx.OrientationBP2 < 0).add(ctx.OrientationBP2 > 0 ? 0 : -1))
                    .filter(l -> withinRange(l, ctx.MantaBP2, ctx.Uncertainty2, extraUncertainty))
                    .max((a, b) -> ctx.OrientationBP2 > 0 ? a.compareTo(b) : b.compareTo(a))
                    .orElse(null);

            if (candidate == null) {
                return BreakpointResult.from(BreakpointError.NO_EVIDENCE);
            }

            bp2_candidates.add(candidate);
        }

        Location breakpoint1 = null;
        Location breakpoint2 = null;
        long max_count = 0;

        for (final Location candidate1 : bp1_candidates) {
            for (final Location candidate2 : bp2_candidates) {

                if (candidate1.compareTo(candidate2) >= 0) {
                    continue;
                }

                final long count = interesting.stream().
                        filter(p -> ctx.OrientationBP1 > 0
                                ? Location.fromSAMRecord(p.getLeft(), true).compareTo(candidate1) <= 0
                                : Location.fromSAMRecord(p.getLeft(), false).compareTo(candidate1) >= 0).
                        filter(p -> ctx.OrientationBP2 > 0
                                ? Location.fromSAMRecord(p.getRight(), true).compareTo(candidate2) <= 0
                                : Location.fromSAMRecord(p.getRight(), false).compareTo(candidate2) >= 0).count();

                if (count > max_count) {
                    breakpoint1 = ctx.OrientationBP1 > 0 ? candidate1.add(-1) : candidate1;
                    breakpoint2 = ctx.OrientationBP2 > 0 ? candidate2.add(-1) : candidate2;
                    max_count = count;
                }

            }
        }

        if (breakpoint1 == null || breakpoint2 == null) {
            return BreakpointResult.from(BreakpointError.ALGO_ERROR);
        }

        return BreakpointResult.from(Pair.of(breakpoint1, breakpoint2));
    }

    private static PairedReads pairs(final AlignmentMap alignments) {
        final PairedReads pairs = new PairedReads(alignments.size());
        for (final AlignmentList list : alignments.values()) {
            for (int i = 0; i < list.size(); ++i) {

                final SAMRecord r0 = list.get(i);
                if (r0.getReadUnmappedFlag()) {
                    continue;
                }

                for (int j = i + 1; j < list.size(); ++j) {

                    final SAMRecord r1 = list.get(j);
                    if (r1.getReadUnmappedFlag()) {
                        continue;
                    }

                    if (isMate(r0, r1) || isMate(r1, r0)) {
                        pairs.add(Pair.of(r0, r1));
                    }

                }

            }
        }
        return pairs;
    }

    private static AlignmentMap readsByName(final List<SAMRecord> alignments) {
        final AlignmentMap result = new AlignmentMap(alignments.size());
        alignments.forEach(a -> result.computeIfAbsent(a.getReadName(), k -> new AlignmentList()).add(a));
        return result;
    }

    private static List<SAMRecord> performQuery(final SamReader reader, final QueryInterval[] intervals) {
        return reader.queryOverlapping(intervals).toList();
    }

    StructuralVariantResult processStructuralVariant(final HMFVariantContext ctx) {

        // perform query for reads

        // TODO: can we query a tighter range?
        QueryInterval[] queryIntervals =
                { new QueryInterval(ctx.MantaBP1.ReferenceIndex, Math.max(0, ctx.MantaBP1.Position + ctx.Uncertainty1.Start - range),
                        ctx.MantaBP1.Position + ctx.Uncertainty1.End + range),
                        new QueryInterval(ctx.MantaBP2.ReferenceIndex, Math.max(0, ctx.MantaBP2.Position + ctx.Uncertainty2.Start - range),
                                ctx.MantaBP2.Position + ctx.Uncertainty2.End + range) };
        queryIntervals = QueryInterval.optimizeIntervals(queryIntervals);

        final List<SAMRecord> tumorReads = performQuery(tumorReader, queryIntervals);
        final List<SAMRecord> refReads = performQuery(refReader, queryIntervals);

        // write evidence

        if (tumorWriter != null) {
            for (final SAMRecord read : tumorReads) {
                final int hash = read.getSAMString().hashCode();
                if (!tumorWrittenReads.contains(hash)) {
                    tumorWriter.addAlignment(read);
                    tumorWrittenReads.add(hash);
                }
            }
        }
        if (refWriter != null) {
            for (final SAMRecord read : refReads) {
                final int hash = read.getSAMString().hashCode();
                if (!refWrittenReads.contains(hash)) {
                    refWriter.addAlignment(read);
                    refWrittenReads.add(hash);
                }
            }
        }

        // process reads

        final PairedReads tumorPairedReads = pairs(readsByName(tumorReads));
        final PairedReads refPairedReads = pairs(readsByName(refReads));

        int extraUncertainty = 0;
        BreakpointResult breakpoints = new BreakpointResult(BreakpointError.UNINITIALIZED);
        for (final int u : extraUncertainties) {
            extraUncertainty = u;
            breakpoints = determineBreakpoints(ctx, tumorPairedReads, u);
            if (breakpoints.Error == BreakpointError.NONE) {
                break;
            }
        }

        final StructuralVariantResult result = new StructuralVariantResult();
        result.Breakpoints = breakpoints.Breakpoints;
        result.ExtraUncertainty = extraUncertainty;

        if (breakpoints.Error != BreakpointError.NONE) {
            LOGGER.error("breakpoint error : {}", ctx.Id);
            result.Filters = Filter.getErrorFilter();
        } else {
            result.TumorStats = collectEvidence(ctx, tumorPairedReads, result.Breakpoints);
            result.RefStats = collectEvidence(ctx, refPairedReads, result.Breakpoints);
            result.AlleleFrequency = AlleleFrequency.calculate(ctx, result.TumorStats);

            // load sample clipping
            tumorReads.forEach(r -> Clipping.getClips(r).forEach(c -> result.TumorStats.Sample_Clipping.add(c)));
            refReads.forEach(r -> Clipping.getClips(r).forEach(c -> result.RefStats.Sample_Clipping.add(c)));

            result.Filters = Filter.getFilters(ctx, result.TumorStats, result.RefStats, result.Breakpoints);
        }

        result.FilterString = result.Filters.isEmpty() ? "PASS" : String.join(";", result.Filters);

        return result;
    }
}
