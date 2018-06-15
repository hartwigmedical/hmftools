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
import com.hartwig.hmftools.breakpointinspector.datamodel.Range;

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

    private static boolean overlap(final Pair<SAMRecord, SAMRecord> pair, final Location breakpoint) {
        return stream(pair).anyMatch(r -> overlap(r, breakpoint));
    }

    private static boolean clippedOnCorrectSide(final SAMRecord record, final int orientation) {
        return (orientation > 0 ? record.getCigar().isRightClipped() : record.getCigar().isLeftClipped());
    }

    private static boolean exactlyClipsBreakpoint(final SAMRecord record, final Location breakpoint, final int orientation) {
        return Location.fromSAMRecord(record, orientation < 0).compareTo(breakpoint) == 0 && clippedOnCorrectSide(record, orientation);
    }

    private boolean withinRange(final Location a, final Location b, final Range range) {
        final int extraUncertainty = 1;
        return a.ReferenceIndex == b.ReferenceIndex && (a.Position >= b.Position + range.start() - extraUncertainty) && (a.Position
                <= b.Position + range.end() + extraUncertainty);
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

        final Location bp1 = breakpoints.getLeft();
        final Location bp2 = breakpoints.getRight();

        final SampleStats result = new SampleStats();
        final Pair<Integer, Integer> ctxOrientation = Pair.of(ctx.OrientationBP1, ctx.OrientationBP2);

        final List<SAMRecord> currentReads = Lists.newArrayList();
        final SAMRecordIterator iterator = reader.iterator();

        final boolean srOnly = ctx.isShortVariant() || ctx.isInsert();

        // iterate through all records in the bam
        // then go through alignments of a read pair-wise
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

            boolean pr_support = false;
            boolean bp1_sr_support = false;
            boolean bp2_sr_support = false;

            boolean bp1_pr_normal = false;
            boolean bp1_sr_normal = false;
            boolean bp2_pr_normal = false;
            boolean bp2_sr_normal = false;

            for (final Pair<SAMRecord, SAMRecord> pair : pairs) {

                final boolean proper = stream(pair).anyMatch(SAMRecord::getProperPairFlag);
                final boolean secondary = stream(pair).anyMatch(SAMRecord::isSecondaryOrSupplementary);

                final boolean correctOrientation = orientation(pair).equals(ctxOrientation);
                final boolean correctChromosome =
                        Location.fromSAMRecord(pair.getLeft()).sameChromosomeAs(ctx.MantaBP1) && Location.fromSAMRecord(pair.getRight())
                                .sameChromosomeAs(ctx.MantaBP2);

                final int MAX_INTRA_PAIR_LENGTH = 400;
                final boolean intraPairLength = (ctx.OrientationBP1 > 0
                        ? bp1.Position - pair.getLeft().getAlignmentEnd()
                        : pair.getLeft().getAlignmentStart() - bp1.Position) + (ctx.OrientationBP2 > 0 ? breakpoints.getRight().Position
                        - pair.getRight().getAlignmentEnd() : pair.getRight().getAlignmentStart() - bp2.Position) < MAX_INTRA_PAIR_LENGTH;

                boolean isPairEvidence = correctOrientation && correctChromosome && intraPairLength;
                if (isPairEvidence) {

                    final int left_outer = Location.fromSAMRecord(pair.getLeft(), ctx.OrientationBP1 > 0).compareTo(bp1);
                    final int right_outer = Location.fromSAMRecord(pair.getRight(), ctx.OrientationBP2 > 0).compareTo(bp2);

                    if (ctx.OrientationBP1 > 0) {
                        isPairEvidence &= left_outer < 0;
                    } else {
                        isPairEvidence &= left_outer > 0;
                    }
                    if (ctx.OrientationBP2 > 0) {
                        isPairEvidence &= right_outer < 0;
                    } else {
                        isPairEvidence &= right_outer > 0;
                    }

                }

                if (isPairEvidence) {
                    bp1_sr_support |= exactlyClipsBreakpoint(pair.getLeft(), bp1, ctx.OrientationBP1);
                    bp2_sr_support |= exactlyClipsBreakpoint(pair.getRight(), bp2, ctx.OrientationBP2);
                    if (!srOnly) {
                        pr_support = true;
                        result.PR_Evidence.add(pair);
                    }
                }

                if (proper || secondary) {

                    final boolean clips_bp1 =
                            exactlyClipsBreakpoint(ctx.OrientationBP1 > 0 ? pair.getRight() : pair.getLeft(), bp1, ctx.OrientationBP1);
                    final boolean clips_bp2 =
                            exactlyClipsBreakpoint(ctx.OrientationBP2 > 0 ? pair.getRight() : pair.getLeft(), bp2, ctx.OrientationBP2);

                    final boolean span_bp1 = span(pair, bp1);
                    final boolean span_bp2 = span(pair, bp2);
                    final boolean overlap_bp1 = overlap(pair, bp1);
                    final boolean overlap_bp2 = overlap(pair, bp2);

                    boolean addToSR = false;
                    if (span_bp1) {
                        if (clips_bp1) {
                            bp1_sr_support = addToSR = true;
                        } else if (overlap_bp1) {
                            bp1_pr_normal = bp1_sr_normal = addToSR = true;
                        } else {
                            bp1_pr_normal = true;
                        }
                    }
                    if (span_bp2) {
                        if (clips_bp2) {
                            bp2_sr_support = addToSR = true;
                        } else if (overlap_bp2) {
                            bp2_pr_normal = bp2_sr_normal = addToSR = true;
                        } else {
                            bp2_pr_normal = true;
                        }
                    }

                    if (addToSR) {
                        result.SR_Evidence.add(pair);
                    }
                }

            } // next pair in reads

            // increment read counts
            final boolean sr_support = bp1_sr_support || bp2_sr_support;

            if (sr_support && pr_support) {
                result.BP1_Stats.PR_SR_Support++;
            } else if (bp1_sr_support) {
                result.BP1_Stats.SR_Only_Support++;
            } else if (pr_support) {
                result.BP1_Stats.PR_Only_Support++;
            }
            if (bp1_pr_normal && bp1_sr_normal) {
                result.BP1_Stats.PR_SR_Normal++;
            } else if (bp1_pr_normal && !srOnly) {
                result.BP1_Stats.PR_Only_Normal++;
            }

            if (sr_support && pr_support) {
                result.BP2_Stats.PR_SR_Support++;
            } else if (bp2_sr_support) {
                result.BP2_Stats.SR_Only_Support++;
            } else if (pr_support) {
                result.BP2_Stats.PR_Only_Support++;
            }
            if (bp2_pr_normal && bp2_sr_normal) {
                result.BP2_Stats.PR_SR_Normal++;
            } else if (bp2_pr_normal && !srOnly) {
                result.BP2_Stats.PR_Only_Normal++;
            }

        } // next read collection

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
            final Location bp2 = ctx.MantaBP2.add(ctx.OrientationBP2 > 0 ? ctx.Uncertainty2.end() : ctx.Uncertainty2.start() + adj);
            return BreakpointResult.from(Pair.of(bp1, bp2));
        } else {
            final Location bp1 = ctx.MantaBP1.add(ctx.OrientationBP1 > 0 ? 0 : adj); // ignore homology when we have an insert
            // TODO: double check adding uncertainty on bp2?
            final Location bp2 = ctx.MantaBP2.add(ctx.OrientationBP2 > 0 ? ctx.Uncertainty2.end() : ctx.Uncertainty2.start() + adj);
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
                new QueryInterval(ctx.MantaBP1.ReferenceIndex, Math.max(0, ctx.MantaBP1.Position + ctx.Uncertainty1.start() - range),
                        ctx.MantaBP1.Position + ctx.Uncertainty1.end() + range),
                new QueryInterval(ctx.MantaBP2.ReferenceIndex, Math.max(0, ctx.MantaBP2.Position + ctx.Uncertainty2.start() - range),
                        ctx.MantaBP2.Position + ctx.Uncertainty2.end() + range) });

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
            result.AlleleFrequency = AlleleFrequency.calculate(result.TumorStats);

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
