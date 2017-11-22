package com.hartwig.hmftools.breakpointinspector;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.breakpointinspector.clipping.Clipping;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

class Analysis {

    private static final Logger LOGGER = LogManager.getLogger(Analysis.class);

    private final SamReader refReader;
    private final SamReader tumorReader;

    private final int range;
    private final float contamination;

    private static final SAMRecordCoordinateComparator comparator = new SAMRecordCoordinateComparator();

    Analysis(final SamReader refReader, final SamReader tumorReader, final int range, final float contamination) {
        this.refReader = refReader;
        this.tumorReader = tumorReader;
        this.range = range;
        this.contamination = contamination;
    }

    private static class PairedReads extends ArrayList<Pair<SAMRecord, SAMRecord>> {
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

    private static boolean exactlyClipsBreakpoint(final SAMRecord record, final Location breakpoint, final int orientation) {
        return Location.fromSAMRecord(record, orientation < 0).compareTo(breakpoint) == 0 && clippedOnCorrectSide(record, orientation);
    }

    private boolean withinRange(final Location a, final Location b, final Range range) {
        final int extraUncertainty = 1;
        return a.ReferenceIndex == b.ReferenceIndex && (a.Position >= b.Position + range.Start - extraUncertainty) && (a.Position
                <= b.Position + range.End + extraUncertainty);
    }

    private static PairedReads pairs(final List<SAMRecord> list) {
        final PairedReads pairs = new PairedReads();
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

                // check both directions due to secondary alignments
                if (isMate(r0, r1) || isMate(r1, r0)) {
                    pairs.add(Pair.of(r0, r1));
                }

            }

        }
        return pairs;
    }

    private SampleStats collectEvidence(final HMFVariantContext ctx, final SamReader reader, final Pair<Location, Location> breakpoints) {

        final SampleStats result = new SampleStats();
        final Pair<Integer, Integer> ctxOrientation = Pair.of(ctx.OrientationBP1, ctx.OrientationBP2);

        final List<SAMRecord> currentReads = Lists.newArrayList();
        final SAMRecordIterator iterator = reader.iterator();

        while (iterator.hasNext() || !currentReads.isEmpty()) {

            final SAMRecord record = iterator.hasNext() ? iterator.next() : null;
            if (record != null) {
                if (currentReads.isEmpty() || record.getReadName().equals(currentReads.get(0).getReadName())) {
                    currentReads.add(record);
                    continue;
                }
            }

            currentReads.sort(comparator);
            final PairedReads pairs = pairs(currentReads);

            currentReads.clear();
            if (record != null) {
                currentReads.add(record);
            }

            for (final Pair<SAMRecord, SAMRecord> pair : pairs) {

                final boolean proper = stream(pair).anyMatch(SAMRecord::getProperPairFlag);
                final boolean secondary = stream(pair).anyMatch(SAMRecord::isSecondaryOrSupplementary);
                final boolean correctOrientation = orientation(pair).equals(ctxOrientation);
                final boolean correctChromosome =
                        Location.fromSAMRecord(pair.getLeft()).sameChromosomeAs(ctx.MantaBP1) && Location.fromSAMRecord(pair.getRight())
                                .sameChromosomeAs(ctx.MantaBP2);

                final int MAX_INTRA_PAIR_LENGTH = 400;
                final boolean intraPairLength = (ctx.OrientationBP1 > 0
                        ? breakpoints.getLeft().Position - pair.getLeft().getAlignmentEnd()
                        : pair.getLeft().getAlignmentStart() - breakpoints.getLeft().Position) + (ctx.OrientationBP2 > 0
                        ? breakpoints.getRight().Position - pair.getRight().getAlignmentEnd()
                        : pair.getRight().getAlignmentStart() - breakpoints.getRight().Position) < MAX_INTRA_PAIR_LENGTH;

                LOGGER.trace("collectEvidence {} {}->{} {} {}->{} correctOrientation({}) correctChromosome({}) correctPosition({})",
                        pair.getLeft().getReadName(), pair.getLeft().getAlignmentStart(), pair.getLeft().getMateAlignmentStart(),
                        pair.getRight().getReadName(), pair.getRight().getAlignmentStart(), pair.getRight().getMateAlignmentStart(),
                        correctOrientation, correctChromosome, intraPairLength);

                boolean support = correctOrientation && correctChromosome && intraPairLength;
                if (support) {

                    final int left_outer = Location.fromSAMRecord(pair.getLeft(), ctx.OrientationBP1 > 0).compareTo(breakpoints.getLeft());
                    final int left_inner = Location.fromSAMRecord(pair.getLeft(), ctx.OrientationBP1 < 0).compareTo(breakpoints.getLeft());
                    final int right_outer =
                            Location.fromSAMRecord(pair.getRight(), ctx.OrientationBP2 > 0).compareTo(breakpoints.getRight());
                    final int right_inner =
                            Location.fromSAMRecord(pair.getRight(), ctx.OrientationBP2 < 0).compareTo(breakpoints.getRight());

                    LOGGER.trace("collectEvidence {} {}->{} {} {}->{} leftOuter({}) leftInner({}) rightInner({}) rightOuter({})",
                            pair.getLeft().getReadName(), pair.getLeft().getAlignmentStart(), pair.getLeft().getMateAlignmentStart(),
                            pair.getRight().getReadName(), pair.getRight().getAlignmentStart(), pair.getRight().getMateAlignmentStart(),
                            left_outer, left_inner, right_inner, right_outer);

                    if (ctx.OrientationBP1 > 0) {
                        support &= left_outer < 0;
                        //support &= left_inner <= 0;
                    } else {
                        support &= left_outer > 0;
                        //support &= left_inner >= 0;
                    }
                    if (ctx.OrientationBP2 > 0) {
                        support &= right_outer < 0;
                        //support &= right_inner <= 0;
                    } else {
                        support &= right_outer > 0;
                        //support &= right_inner >= 0;
                    }

                }

                if (support) {

                    final boolean clip_bp1 = exactlyClipsBreakpoint(pair.getLeft(), breakpoints.getLeft(), ctx.OrientationBP1);
                    final boolean clip_bp2 = exactlyClipsBreakpoint(pair.getRight(), breakpoints.getRight(), ctx.OrientationBP2);

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

                }

                if (proper || secondary) {

                    final boolean clip_bp1 =
                            exactlyClipsBreakpoint(ctx.OrientationBP1 > 0 ? pair.getRight() : pair.getLeft(), breakpoints.getLeft(),
                                    ctx.OrientationBP1);
                    final boolean clip_bp2 =
                            exactlyClipsBreakpoint(ctx.OrientationBP2 > 0 ? pair.getRight() : pair.getLeft(), breakpoints.getRight(),
                                    ctx.OrientationBP2);

                    final boolean span_bp1 = span(pair, breakpoints.getLeft());
                    final boolean span_bp2 = span(pair, breakpoints.getRight());
                    final boolean overlap_bp1 = stream(pair).anyMatch(r -> overlap(r, breakpoints.getLeft()));
                    final boolean overlap_bp2 = stream(pair).anyMatch(r -> overlap(r, breakpoints.getRight()));

                    if (clip_bp1) {
                        result.BP1_Stats.SR_Only_Support++;
                        result.SR_Evidence.add(pair);
                        continue;
                    }
                    if (clip_bp2) {
                        result.BP2_Stats.SR_Only_Support++;
                        result.SR_Evidence.add(pair);
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
        }

        iterator.close();
        return result;
    }

    enum BreakpointError {
        NONE,
        ALGO_ERROR
    }

    private static class BreakpointResult {
        private BreakpointResult(final Pair<Location, Location> breakpoints) {
            Breakpoints = breakpoints;
            if (stream(breakpoints).anyMatch(Objects::isNull)) {
                Error = BreakpointError.ALGO_ERROR;
            }
        }

        static BreakpointResult from(final Pair<Location, Location> breakpoints) {
            return new BreakpointResult(breakpoints);
        }

        Pair<Location, Location> Breakpoints;
        BreakpointError Error = BreakpointError.NONE;
    }

    private BreakpointResult determineBreakpointsImprecise(final HMFVariantContext ctx, final SamReader reader) {

        final Pair<Integer, Integer> ctxOrientation = Pair.of(ctx.OrientationBP1, ctx.OrientationBP2);

        final PairedReads interesting = new PairedReads();
        final PairedReads clipped_proper = new PairedReads();
        final PairedReads secondary_pairs = new PairedReads();

        final List<SAMRecord> currentReads = Lists.newArrayList();
        final SAMRecordIterator iterator = reader.iterator();

        while (iterator.hasNext() || !currentReads.isEmpty()) {

            final SAMRecord record = iterator.hasNext() ? iterator.next() : null;
            if (record != null) {
                if (currentReads.isEmpty() || record.getReadName().equals(currentReads.get(0).getReadName())) {
                    currentReads.add(record);
                    continue;
                }
            }

            currentReads.sort(comparator);
            final PairedReads pairs = pairs(currentReads);

            currentReads.clear();
            if (record != null) {
                currentReads.add(record);
            }

            // extract all interesting pairs

            for (final Pair<SAMRecord, SAMRecord> pair : pairs) {

                final boolean correctOrientation = orientation(pair).equals(ctxOrientation);
                final boolean correctChromosome =
                        Location.fromSAMRecord(pair.getLeft()).sameChromosomeAs(ctx.MantaBP1) && Location.fromSAMRecord(pair.getRight())
                                .sameChromosomeAs(ctx.MantaBP2);
                final boolean hasExpectedClipping =
                        clippedOnCorrectSide(pair.getLeft(), ctx.OrientationBP1) || clippedOnCorrectSide(pair.getRight(),
                                ctx.OrientationBP2);

                final boolean sameChromosome = pair.getLeft().getReferenceIndex().equals(pair.getRight().getReferenceIndex());
                final boolean potentialSROnly = sameChromosome && Stream.of(ctx.OrientationBP1, ctx.OrientationBP2)
                        .anyMatch(orientation -> clippedOnCorrectSide(orientation > 0 ? pair.getRight() : pair.getLeft(), orientation));

                final boolean secondary = stream(pair).anyMatch(SAMRecord::isSecondaryOrSupplementary);
                final boolean proper = stream(pair).anyMatch(SAMRecord::getProperPairFlag) && !secondary;

                LOGGER.trace(
                        "determineBreakpoints {} {}->{} {} {}->{} correctOrientation({}) correctChromosome({}) hasExpectedClipping({}) proper({}) potentialSROnly({}) secondary({})",
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

        }

        iterator.close();

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

        final List<Location> bp1_candidates = bp1_clipping.getSequences().stream()
                //.filter(c -> c.LongestClipSequence.length() >= 5)
                .map(c -> c.Alignment).filter(c -> withinRange(c, ctx.MantaBP1, ctx.Uncertainty1)).collect(Collectors.toList());

        if (bp1_candidates.isEmpty()) {
            final Location candidate = interesting.stream()
                    .map(Pair::getLeft)
                    .map(r -> Location.fromSAMRecord(r, ctx.OrientationBP1 < 0).add(ctx.OrientationBP1 > 0 ? 1 : -1))
                    .filter(l -> withinRange(l, ctx.MantaBP1, ctx.Uncertainty1))
                    .max((a, b) -> ctx.OrientationBP1 > 0 ? a.compareTo(b) : b.compareTo(a))
                    .orElse(null);

            if (candidate != null) {
                bp1_candidates.add(candidate);
            }
        }

        final List<Location> bp2_candidates = bp2_clipping.getSequences().stream()
                //.filter(c -> c.LongestClipSequence.length() >= 5)
                .map(c -> c.Alignment).filter(c -> withinRange(c, ctx.MantaBP2, ctx.Uncertainty2)).collect(Collectors.toList());

        if (bp2_candidates.isEmpty()) {
            final Location candidate = interesting.stream()
                    .map(Pair::getRight)
                    .map(r -> Location.fromSAMRecord(r, ctx.OrientationBP2 < 0).add(ctx.OrientationBP2 > 0 ? 1 : -1))
                    .filter(l -> withinRange(l, ctx.MantaBP2, ctx.Uncertainty2))
                    .max((a, b) -> ctx.OrientationBP2 > 0 ? a.compareTo(b) : b.compareTo(a))
                    .orElse(null);

            if (candidate != null) {
                bp2_candidates.add(candidate);
            }
        }

        // NOTE: we include homology on both sides here and take it out later
        LOGGER.trace("bp1_candidates={} bp2_candidates={}", bp1_candidates, bp2_candidates);
        final Location breakpoint1 = bp1_candidates.isEmpty() ? null : bp1_candidates.get(0).add(-ctx.OrientationBP1);
        final Location breakpoint2 = bp2_candidates.isEmpty() ? null : bp2_candidates.get(0).add(-ctx.OrientationBP2);

        return BreakpointResult.from(Pair.of(breakpoint1, breakpoint2));
    }

    // basically, this will align to where we expect to see clipping
    private BreakpointResult determineBreakpoints(final HMFVariantContext ctx, final SamReader reader) {
        final int adj = ctx.BND ? 0 : 1;
        if (ctx.Imprecise) {
            return determineBreakpointsImprecise(ctx, reader);
        } else if (ctx.isInsert()) {
            return BreakpointResult.from(Pair.of(ctx.MantaBP1, ctx.MantaBP2.add(1))); // we want last match base at this stage
        } else if (ctx.InsertSequence.isEmpty()) {
            final Location bp1 = ctx.MantaBP1.add(ctx.OrientationBP1 > 0 ? ctx.HomologySequence.length() : adj);
            final Location bp2 = ctx.MantaBP2.add(ctx.OrientationBP2 > 0 ? ctx.Uncertainty2.End : ctx.Uncertainty2.Start + adj);
            return BreakpointResult.from(Pair.of(bp1, bp2));
        } else {
            final Location bp1 = ctx.MantaBP1.add(ctx.OrientationBP1 > 0 ? 0 : adj); // ignore homology when we have an insert
            final Location bp2 = ctx.MantaBP2.add(ctx.OrientationBP2 > 0 ? ctx.Uncertainty2.End : ctx.Uncertainty2.Start + adj);
            return BreakpointResult.from(Pair.of(bp1, bp2));
        }
    }

    private static File queryNameSortedBAM(final SamReader reader, final QueryInterval[] intervals, final String name) throws IOException {

        final SAMFileHeader header = reader.getFileHeader().clone();
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);

        final File file = File.createTempFile(name, ".bam");
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, file);

        final SAMRecordIterator iterator = reader.queryOverlapping(intervals);
        while (iterator.hasNext()) {
            writer.addAlignment(iterator.next());
        }

        iterator.close();
        writer.close();

        return file;
    }

    StructuralVariantResult processStructuralVariant(final HMFVariantContext ctx) throws IOException {

        final QueryInterval[] intervals = QueryInterval.optimizeIntervals(new QueryInterval[] {
                new QueryInterval(ctx.MantaBP1.ReferenceIndex, Math.max(0, ctx.MantaBP1.Position + ctx.Uncertainty1.Start - range),
                        ctx.MantaBP1.Position + ctx.Uncertainty1.End + range),
                new QueryInterval(ctx.MantaBP2.ReferenceIndex, Math.max(0, ctx.MantaBP2.Position + ctx.Uncertainty2.Start - range),
                        ctx.MantaBP2.Position + ctx.Uncertainty2.End + range) });

        final File TEMP_REF_BAM = queryNameSortedBAM(refReader, intervals, "ref");
        final File TEMP_TUMOR_BAM = queryNameSortedBAM(tumorReader, intervals, "tumor");

        final SamReader SORTED_REF_READER = SamReaderFactory.makeDefault().open(TEMP_REF_BAM);
        final SamReader SORTED_TUMOR_READER = SamReaderFactory.makeDefault().open(TEMP_TUMOR_BAM);

        final BreakpointResult breakpoints = determineBreakpoints(ctx, SORTED_TUMOR_READER);

        final StructuralVariantResult result = new StructuralVariantResult();
        result.Breakpoints = breakpoints.Breakpoints;
        result.QueryIntervals = intervals;

        if (breakpoints.Error != BreakpointError.NONE) {
            result.Filters = Filter.getErrorFilter();
        } else {
            result.TumorStats = collectEvidence(ctx, SORTED_TUMOR_READER, result.Breakpoints);
            result.RefStats = collectEvidence(ctx, SORTED_REF_READER, result.Breakpoints);
            result.AlleleFrequency = AlleleFrequency.calculate(ctx, result.TumorStats);

            // load sample clipping
            SORTED_TUMOR_READER.forEach(r -> Clipping.getClips(r).forEach(c -> result.TumorStats.Sample_Clipping.add(c)));
            SORTED_REF_READER.forEach(r -> Clipping.getClips(r).forEach(c -> result.RefStats.Sample_Clipping.add(c)));

            result.Filters = Filter.getFilters(ctx, result.TumorStats, result.RefStats, result.Breakpoints, contamination);

            // adjust for homology
            final Location bp1 = result.Breakpoints.getLeft().add(ctx.OrientationBP1 > 0 ? 0 : -1);
            final Location bp2;
            if (!ctx.isInsert() && ctx.InsertSequence.isEmpty()) {
                bp2 = result.Breakpoints.getRight()
                        .add(-ctx.OrientationBP2 * ctx.HomologySequence.length())
                        .add(ctx.OrientationBP2 > 0 ? 0 : -1);
            } else {
                bp2 = result.Breakpoints.getRight().add(ctx.OrientationBP2 > 0 ? 0 : -1);
            }
            result.Breakpoints = Pair.of(bp1, bp2);
        }

        result.FilterString = result.Filters.isEmpty() ? "PASS" : String.join(";", result.Filters);

        // clean up
        SORTED_REF_READER.close();
        SORTED_TUMOR_READER.close();

        if (!TEMP_REF_BAM.delete()) {
            LOGGER.error("couldn't delete {}", TEMP_REF_BAM);
        }
        if (!TEMP_TUMOR_BAM.delete()) {
            LOGGER.error("couldn't delete {}", TEMP_TUMOR_BAM);
        }

        return result;
    }
}
