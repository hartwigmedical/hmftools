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
import com.hartwig.hmftools.breakpointinspector.datamodel.EnrichedVariantContext;
import com.hartwig.hmftools.breakpointinspector.datamodel.Range;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

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
    private static final SAMRecordCoordinateComparator COMPARATOR = new SAMRecordCoordinateComparator();

    private final SamReader refReader;
    private final SamReader tumorReader;

    private final int range;
    private final float contamination;

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
        return read.getReferenceIndex() == breakpoint.referenceIndex() && read.getAlignmentStart() <= breakpoint.position()
                && breakpoint.position() <= read.getAlignmentEnd();
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

    private static boolean withinRange(final Location a, final Location b, final Range range) {
        final int extraUncertainty = 1;
        return a.referenceIndex() == b.referenceIndex() && (a.position() >= b.position() + range.start() - extraUncertainty) && (
                a.position() <= b.position() + range.end() + extraUncertainty);
    }

    @NotNull
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

    @NotNull
    private static SampleStats collectEvidence(final EnrichedVariantContext context, final SamReader reader,
            final Pair<Location, Location> breakpoints) {
        final Location bp1 = breakpoints.getLeft();
        final Location bp2 = breakpoints.getRight();

        final SampleStats result = new SampleStats();
        final Pair<Integer, Integer> ctxOrientation = Pair.of(context.orientationBP1(), context.orientationBP2());

        final List<SAMRecord> currentReads = Lists.newArrayList();
        final SAMRecordIterator iterator = reader.iterator();

        final boolean srOnly = context.isShortVariant() || context.isInsert();

        // NERA: Iterate through all records in the bam, then go through alignments of a read pair-wise
        while (iterator.hasNext() || !currentReads.isEmpty()) {
            final SAMRecord record = iterator.hasNext() ? iterator.next() : null;
            if (record != null) {
                if (currentReads.isEmpty() || record.getReadName().equals(currentReads.get(0).getReadName())) {
                    currentReads.add(record);
                    continue;
                }
            }

            currentReads.sort(COMPARATOR);
            final PairedReads pairs = pairs(currentReads);

            currentReads.clear();
            if (record != null) {
                currentReads.add(record);
            }

            boolean prSupport = false;
            boolean bp1SRSupport = false;
            boolean bp2SRSupport = false;

            boolean bp1PRNormal = false;
            boolean bp1SRNormal = false;
            boolean bp2PRNormal = false;
            boolean bp2SRNormal = false;

            for (final Pair<SAMRecord, SAMRecord> pair : pairs) {
                final boolean proper = stream(pair).anyMatch(SAMRecord::getProperPairFlag);
                final boolean secondary = stream(pair).anyMatch(SAMRecord::isSecondaryOrSupplementary);

                final boolean correctOrientation = orientation(pair).equals(ctxOrientation);
                final boolean correctChromosome = Location.fromSAMRecord(pair.getLeft()).sameChromosomeAs(context.locationBP1()) && Location
                        .fromSAMRecord(pair.getRight())
                        .sameChromosomeAs(context.locationBP2());

                final int MAX_INTRA_PAIR_LENGTH = 400;
                final boolean intraPairLength = (context.orientationBP1() > 0
                        ? bp1.position() - pair.getLeft().getAlignmentEnd()
                        : pair.getLeft().getAlignmentStart() - bp1.position()) + (context.orientationBP2() > 0
                        ? breakpoints.getRight().position() - pair.getRight().getAlignmentEnd()
                        : pair.getRight().getAlignmentStart() - bp2.position()) < MAX_INTRA_PAIR_LENGTH;

                boolean isPairEvidence = correctOrientation && correctChromosome && intraPairLength;
                if (isPairEvidence) {
                    final int leftOuter = Location.fromSAMRecord(pair.getLeft(), context.orientationBP1() > 0).compareTo(bp1);
                    final int rightOuter = Location.fromSAMRecord(pair.getRight(), context.orientationBP2() > 0).compareTo(bp2);

                    if (context.orientationBP1() > 0) {
                        isPairEvidence &= leftOuter < 0;
                    } else {
                        isPairEvidence &= leftOuter > 0;
                    }
                    if (context.orientationBP2() > 0) {
                        isPairEvidence &= rightOuter < 0;
                    } else {
                        isPairEvidence &= rightOuter > 0;
                    }
                }

                if (isPairEvidence) {
                    bp1SRSupport |= exactlyClipsBreakpoint(pair.getLeft(), bp1, context.orientationBP1());
                    bp2SRSupport |= exactlyClipsBreakpoint(pair.getRight(), bp2, context.orientationBP2());
                    if (!srOnly) {
                        prSupport = true;
                        result.PR_Evidence.add(pair);
                    }
                }

                if (proper || secondary) {
                    final boolean clipsBP1 = exactlyClipsBreakpoint(context.orientationBP1() > 0 ? pair.getRight() : pair.getLeft(),
                            bp1,
                            context.orientationBP1());
                    final boolean clipsBP2 = exactlyClipsBreakpoint(context.orientationBP2() > 0 ? pair.getRight() : pair.getLeft(),
                            bp2,
                            context.orientationBP2());

                    final boolean spanBP1 = span(pair, bp1);
                    final boolean spanBP2 = span(pair, bp2);
                    final boolean overlapBP1 = overlap(pair, bp1);
                    final boolean overlapBP2 = overlap(pair, bp2);

                    boolean addToSR = false;
                    if (spanBP1) {
                        if (clipsBP1) {
                            bp1SRSupport = addToSR = true;
                        } else {
                            bp1PRNormal = !overlapBP1 || (bp1SRNormal = addToSR = true);
                        }
                    }
                    if (spanBP2) {
                        if (clipsBP2) {
                            bp2SRSupport = addToSR = true;
                        } else {
                            bp2PRNormal = !overlapBP2 || (bp2SRNormal = addToSR = true);
                        }
                    }

                    if (addToSR) {
                        result.SR_Evidence.add(pair);
                    }
                }
            }

            // NERA: Increment read counts
            final boolean srSupport = bp1SRSupport || bp2SRSupport;

            if (srSupport && prSupport) {
                result.BP1_Stats.PR_SR_Support++;
            } else if (bp1SRSupport) {
                result.BP1_Stats.SR_Only_Support++;
            } else if (prSupport) {
                result.BP1_Stats.PR_Only_Support++;
            }
            if (bp1PRNormal && bp1SRNormal) {
                result.BP1_Stats.PR_SR_Normal++;
            } else if (bp1PRNormal && !srOnly) {
                result.BP1_Stats.PR_Only_Normal++;
            }

            if (srSupport && prSupport) {
                result.BP2_Stats.PR_SR_Support++;
            } else if (bp2SRSupport) {
                result.BP2_Stats.SR_Only_Support++;
            } else if (prSupport) {
                result.BP2_Stats.PR_Only_Support++;
            }
            if (bp2PRNormal && bp2SRNormal) {
                result.BP2_Stats.PR_SR_Normal++;
            } else if (bp2PRNormal && !srOnly) {
                result.BP2_Stats.PR_Only_Normal++;
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

    private static BreakpointResult determineBreakpointsImprecise(final EnrichedVariantContext variant, final SamReader reader) {
        final Pair<Integer, Integer> ctxOrientation = Pair.of(variant.orientationBP1(), variant.orientationBP2());

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

            currentReads.sort(COMPARATOR);
            final PairedReads pairs = pairs(currentReads);

            currentReads.clear();
            if (record != null) {
                currentReads.add(record);
            }

            // extract all interesting pairs

            for (final Pair<SAMRecord, SAMRecord> pair : pairs) {

                final boolean correctOrientation = orientation(pair).equals(ctxOrientation);
                final boolean correctChromosome = Location.fromSAMRecord(pair.getLeft()).sameChromosomeAs(variant.locationBP1()) && Location
                        .fromSAMRecord(pair.getRight())
                        .sameChromosomeAs(variant.locationBP2());
                final boolean hasExpectedClipping = clippedOnCorrectSide(pair.getLeft(), variant.orientationBP1()) || clippedOnCorrectSide(
                        pair.getRight(),
                        variant.orientationBP2());

                final boolean sameChromosome = pair.getLeft().getReferenceIndex().equals(pair.getRight().getReferenceIndex());
                final boolean potentialSROnly = sameChromosome && Stream.of(variant.orientationBP1(), variant.orientationBP2())
                        .anyMatch(orientation -> clippedOnCorrectSide(orientation > 0 ? pair.getRight() : pair.getLeft(), orientation));

                final boolean secondary = stream(pair).anyMatch(SAMRecord::isSecondaryOrSupplementary);
                final boolean proper = stream(pair).anyMatch(SAMRecord::getProperPairFlag) && !secondary;

                LOGGER.trace(
                        "determineBreakpoints {} {}->{} {} {}->{} correctOrientation({}) correctChromosome({}) hasExpectedClipping({}) proper({}) potentialSROnly({}) secondary({})",
                        pair.getLeft().getReadName(),
                        pair.getLeft().getAlignmentStart(),
                        pair.getLeft().getMateAlignmentStart(),
                        pair.getRight().getReadName(),
                        pair.getRight().getAlignmentStart(),
                        pair.getRight().getMateAlignmentStart(),
                        correctOrientation,
                        correctChromosome,
                        hasExpectedClipping,
                        proper,
                        potentialSROnly,
                        secondary);

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
            if (variant.orientationBP1() > 0) {
                bp1_clipping.add(Clipping.getRightClip(pair.getLeft()));
            } else {
                bp1_clipping.add(Clipping.getLeftClip(pair.getLeft()));
            }
            if (variant.orientationBP2() > 0) {
                bp2_clipping.add(Clipping.getRightClip(pair.getRight()));
            } else {
                bp2_clipping.add(Clipping.getLeftClip(pair.getRight()));
            }
        }

        // include more clipping information

        for (final Pair<SAMRecord, SAMRecord> pair : clipped_proper) {
            if (stream(pair).allMatch(r -> Location.fromSAMRecord(r).sameChromosomeAs(variant.locationBP1()))) {
                if (variant.orientationBP1() > 0) {
                    bp1_clipping.add(Clipping.getRightClip(pair.getRight()));
                } else {
                    bp1_clipping.add(Clipping.getLeftClip(pair.getLeft()));
                }
            }
            if (stream(pair).allMatch(r -> Location.fromSAMRecord(r).sameChromosomeAs(variant.locationBP2()))) {
                if (variant.orientationBP2() > 0) {
                    bp2_clipping.add(Clipping.getRightClip(pair.getRight()));
                } else {
                    bp2_clipping.add(Clipping.getLeftClip(pair.getLeft()));
                }
            }
        }

        // include secondary clipping information

        for (final Pair<SAMRecord, SAMRecord> pair : secondary_pairs) {
            if (stream(pair).allMatch(r -> Location.fromSAMRecord(r).sameChromosomeAs(variant.locationBP1()))) {
                if (variant.orientationBP1() > 0) {
                    bp1_clipping.add(Clipping.getRightClip(pair.getRight()));
                } else {
                    bp1_clipping.add(Clipping.getLeftClip(pair.getLeft()));
                }
            }
            if (stream(pair).allMatch(r -> Location.fromSAMRecord(r).sameChromosomeAs(variant.locationBP2()))) {
                if (variant.orientationBP2() > 0) {
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
                .filter(c -> withinRange(c, variant.locationBP1(), variant.uncertaintyBP1()))
                .collect(Collectors.toList());

        if (bp1_candidates.isEmpty()) {
            interesting.stream()
                    .map(Pair::getLeft)
                    .map(r -> Location.fromSAMRecord(r, variant.orientationBP1() < 0).add(variant.orientationBP1() > 0 ? 1 : -1))
                    .filter(l -> withinRange(l, variant.locationBP1(), variant.uncertaintyBP1()))
                    .max((a, b) -> variant.orientationBP1() > 0 ? a.compareTo(b) : b.compareTo(a))
                    .ifPresent(bp1_candidates::add);
        }

        final List<Location> bp2_candidates = bp2_clipping.getSequences()
                .stream()
                //.filter(c -> c.LongestClipSequence.length() >= 5)
                .map(c -> c.Alignment)
                .filter(c -> withinRange(c, variant.locationBP2(), variant.uncertaintyBP2()))
                .collect(Collectors.toList());

        if (bp2_candidates.isEmpty()) {
            interesting.stream()
                    .map(Pair::getRight)
                    .map(r -> Location.fromSAMRecord(r, variant.orientationBP2() < 0).add(variant.orientationBP2() > 0 ? 1 : -1))
                    .filter(l -> withinRange(l, variant.locationBP2(), variant.uncertaintyBP2()))
                    .max((a, b) -> variant.orientationBP2() > 0 ? a.compareTo(b) : b.compareTo(a))
                    .ifPresent(bp2_candidates::add);
        }

        // NOTE: we include homology on both sides here and take it out later
        LOGGER.trace("bp1_candidates={} bp2_candidates={}", bp1_candidates, bp2_candidates);
        final Location breakpoint1 = bp1_candidates.isEmpty() ? null : bp1_candidates.get(0).add(-variant.orientationBP1());
        final Location breakpoint2 = bp2_candidates.isEmpty() ? null : bp2_candidates.get(0).add(-variant.orientationBP2());

        return BreakpointResult.from(Pair.of(breakpoint1, breakpoint2));
    }

    // basically, this will align to where we expect to see clipping
    private static BreakpointResult determineBreakpoints(final EnrichedVariantContext context, final SamReader reader) {
        final int adj = context.isTranslocation() ? 0 : 1;
        if (context.isImprecise()) {
            return determineBreakpointsImprecise(context, reader);
        } else if (context.isInsert()) {
            return BreakpointResult.from(Pair.of(context.locationBP1(),
                    context.locationBP2().add(1))); // we want last match base at this stage
        } else if (context.insertSequence().isEmpty()) {
            final Location bp1 = context.locationBP1().add(context.orientationBP1() > 0 ? context.homologySequence().length() : adj);
            final Location bp2 = context.locationBP2()
                    .add(context.orientationBP2() > 0 ? context.uncertaintyBP2().end() : context.uncertaintyBP2().start() + adj);
            return BreakpointResult.from(Pair.of(bp1, bp2));
        } else {
            final Location bp1 =
                    context.locationBP1().add(context.orientationBP1() > 0 ? 0 : adj); // ignore homology when we have an insert
            // TODO: double check adding uncertainty on bp2?
            final Location bp2 = context.locationBP2()
                    .add(context.orientationBP2() > 0 ? context.uncertaintyBP2().end() : context.uncertaintyBP2().start() + adj);
            return BreakpointResult.from(Pair.of(bp1, bp2));
        }
    }

    @NotNull
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

    @NotNull
    StructuralVariantResult processStructuralVariant(final EnrichedVariantContext context) throws IOException {
        final QueryInterval[] intervals = QueryInterval.optimizeIntervals(new QueryInterval[] {
                new QueryInterval(context.locationBP1().referenceIndex(),
                        Math.max(0, context.locationBP1().position() + context.uncertaintyBP1().start() - range),
                        context.locationBP1().position() + context.uncertaintyBP1().end() + range),
                new QueryInterval(context.locationBP2().referenceIndex(),
                        Math.max(0, context.locationBP2().position() + context.uncertaintyBP2().start() - range),
                        context.locationBP2().position() + context.uncertaintyBP2().end() + range) });

        final File TEMP_REF_BAM = queryNameSortedBAM(refReader, intervals, "ref");
        final File TEMP_TUMOR_BAM = queryNameSortedBAM(tumorReader, intervals, "tumor");

        final SamReader SORTED_REF_READER = SamReaderFactory.makeDefault().open(TEMP_REF_BAM);
        final SamReader SORTED_TUMOR_READER = SamReaderFactory.makeDefault().open(TEMP_TUMOR_BAM);

        final BreakpointResult breakpoints = determineBreakpoints(context, SORTED_TUMOR_READER);

        final StructuralVariantResult result = new StructuralVariantResult();
        result.Breakpoints = breakpoints.Breakpoints;
        result.QueryIntervals = intervals;

        if (breakpoints.Error != BreakpointError.NONE) {
            result.Filters = Filter.errorFilter();
        } else {
            result.TumorStats = collectEvidence(context, SORTED_TUMOR_READER, result.Breakpoints);
            result.RefStats = collectEvidence(context, SORTED_REF_READER, result.Breakpoints);
            result.AlleleFrequency = AlleleFrequency.calculate(result.TumorStats);

            // NERA: load sample clipping
            SORTED_TUMOR_READER.forEach(r -> Clipping.getClips(r).forEach(c -> result.TumorStats.Sample_Clipping.add(c)));
            SORTED_REF_READER.forEach(r -> Clipping.getClips(r).forEach(c -> result.RefStats.Sample_Clipping.add(c)));

            result.Filters = Filter.filters(context, result.TumorStats, result.RefStats, result.Breakpoints, contamination);

            // NERA: adjust for homology
            final Location bp1 = result.Breakpoints.getLeft().add(context.orientationBP1() > 0 ? 0 : -1);
            final Location bp2;
            if (!context.isInsert() && context.insertSequence().isEmpty()) {
                bp2 = result.Breakpoints.getRight()
                        .add(-context.orientationBP2() * context.homologySequence().length())
                        .add(context.orientationBP2() > 0 ? 0 : -1);
            } else {
                bp2 = result.Breakpoints.getRight().add(context.orientationBP2() > 0 ? 0 : -1);
            }
            result.Breakpoints = Pair.of(bp1, bp2);
        }

        result.FilterString = result.Filters.isEmpty() ? "PASS" : String.join(";", result.Filters);

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
