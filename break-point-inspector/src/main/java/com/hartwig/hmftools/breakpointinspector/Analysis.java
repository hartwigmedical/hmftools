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
import com.hartwig.hmftools.breakpointinspector.datamodel.SampleStats;
import com.hartwig.hmftools.breakpointinspector.datamodel.StructuralVariantResult;

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
    private static final int MAX_READS_PER_INTERVAL_LENGTH_FOR_DOWNSAMPLING_CUTOFF = 100;

    @NotNull
    private final SamReader refReader;
    @NotNull
    private final SamReader tumorReader;

    private final int range;
    private final float contamination;

    Analysis(@NotNull final SamReader refReader, @NotNull final SamReader tumorReader, final int range, final float contamination) {
        this.refReader = refReader;
        this.tumorReader = tumorReader;
        this.range = range;
        this.contamination = contamination;
    }

    @NotNull
    StructuralVariantResult processStructuralVariant(final EnrichedVariantContext variant) throws IOException {
        final QueryInterval[] intervals = QueryInterval.optimizeIntervals(new QueryInterval[] {
                new QueryInterval(variant.locationBP1().referenceIndex(),
                        Math.max(0, variant.locationBP1().position() + variant.uncertaintyBP1().start() - range),
                        variant.locationBP1().position() + variant.uncertaintyBP1().end() + range),
                new QueryInterval(variant.locationBP2().referenceIndex(),
                        Math.max(0, variant.locationBP2().position() + variant.uncertaintyBP2().start() - range),
                        variant.locationBP2().position() + variant.uncertaintyBP2().end() + range) });

        LOGGER.info("  Creating tmp ref BAM using " + Lists.newArrayList(intervals));
        final File tmpRefBam = queryNameSortedBAM(refReader, intervals, "ref");
        LOGGER.info("  Creating tmp tumor BAM using " + Lists.newArrayList(intervals));
        final File tmpTumorBam = queryNameSortedBAM(tumorReader, intervals, "tumor");

        final SamReader sortedRefReader = SamReaderFactory.makeDefault().open(tmpRefBam);
        final SamReader sortedTumorReader = SamReaderFactory.makeDefault().open(tmpTumorBam);

        LOGGER.info("  Analyzing breakpoints");
        final BreakpointResult breakpoints = determineBreakpoints(variant, sortedTumorReader);

        final StructuralVariantResult result = new StructuralVariantResult();
        result.breakpoints = breakpoints.breakpoints;
        result.queryIntervals = intervals;

        if (breakpoints.error != BreakpointError.NONE) {
            result.filters = Filter.errorFilter();
        } else {
            result.tumorStats = collectEvidence(variant, sortedTumorReader, result.breakpoints);
            result.refStats = collectEvidence(variant, sortedRefReader, result.breakpoints);
            result.alleleFrequency = AlleleFrequency.calculate(result.tumorStats);

            // NERA: load sample clipping
            sortedRefReader.forEach(record -> Clipping.clips(record).forEach(clipInfo -> result.refStats.sampleClipping.add(clipInfo)));
            sortedTumorReader.forEach(record -> Clipping.clips(record)
                    .forEach(clipInfo -> result.tumorStats.sampleClipping.add(clipInfo)));

            result.filters = Filter.filters(variant, result.tumorStats, result.refStats, result.breakpoints, contamination);

            // NERA: adjust for homology
            final Location bp1 = result.breakpoints.getLeft().add(variant.orientationBP1() > 0 ? 0 : -1);
            final Location bp2;
            if (!variant.isInsert() && variant.insertSequence().isEmpty()) {
                bp2 = result.breakpoints.getRight()
                        .add(-variant.orientationBP2() * variant.homologySequence().length())
                        .add(variant.orientationBP2() > 0 ? 0 : -1);
            } else {
                bp2 = result.breakpoints.getRight().add(variant.orientationBP2() > 0 ? 0 : -1);
            }
            result.breakpoints = Pair.of(bp1, bp2);
        }

        result.filterString = result.filters.isEmpty() ? "PASS" : String.join(";", result.filters);

        sortedRefReader.close();
        sortedTumorReader.close();

        if (!tmpRefBam.delete()) {
            LOGGER.error(String.format("couldn't delete %s", tmpRefBam));
        }

        if (!tmpTumorBam.delete()) {
            LOGGER.error(String.format("couldn't delete %s", tmpTumorBam));
        }

        return result;
    }

    @NotNull
    private static BreakpointResult determineBreakpoints(final EnrichedVariantContext variant, final SamReader reader) {
        final int adj = variant.isTranslocation() ? 0 : 1;
        if (variant.isImprecise()) {
            return determineBreakpointsImprecise(variant, reader);
        } else if (variant.isInsert()) {
            // NERA: We want last match base at this stage
            return BreakpointResult.from(Pair.of(variant.locationBP1(), variant.locationBP2().add(1)));
        } else if (variant.insertSequence().isEmpty()) {
            final Location bp1 = variant.locationBP1().add(variant.orientationBP1() > 0 ? variant.homologySequence().length() : adj);
            final Location bp2 = variant.locationBP2()
                    .add(variant.orientationBP2() > 0 ? variant.uncertaintyBP2().end() : variant.uncertaintyBP2().start() + adj);
            return BreakpointResult.from(Pair.of(bp1, bp2));
        } else {
            final Location bp1 =
                    variant.locationBP1().add(variant.orientationBP1() > 0 ? 0 : adj); // ignore homology when we have an insert
            final Location bp2 = variant.locationBP2()
                    .add(variant.orientationBP2() > 0 ? variant.uncertaintyBP2().end() : variant.uncertaintyBP2().start() + adj);
            return BreakpointResult.from(Pair.of(bp1, bp2));
        }
    }

    @NotNull
    private static BreakpointResult determineBreakpointsImprecise(@NotNull final EnrichedVariantContext variant,
            @NotNull final SamReader reader) {
        final Pair<Integer, Integer> variantOrientation = Pair.of(variant.orientationBP1(), variant.orientationBP2());

        final PairedReads interesting = new PairedReads();
        final PairedReads clippedProper = new PairedReads();
        final PairedReads secondaryPairs = new PairedReads();

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

            // NERA: extract all interesting pairs
            for (final Pair<SAMRecord, SAMRecord> pair : pairs) {
                final boolean correctOrientation = orientation(pair).equals(variantOrientation);
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

                if (secondary && potentialSROnly) {
                    secondaryPairs.add(pair);
                } else if ((!proper || hasExpectedClipping) && correctChromosome && correctOrientation) {
                    interesting.add(pair);
                } else if (proper && potentialSROnly) {
                    clippedProper.add(pair);
                }
            }

        }

        iterator.close();

        // NERA: load clipping info
        Clipping bp1Clipping = new Clipping();
        Clipping bp2Clipping = new Clipping();

        for (final Pair<SAMRecord, SAMRecord> pair : interesting) {
            if (variant.orientationBP1() > 0) {
                bp1Clipping.add(Clipping.rightClip(pair.getLeft()));
            } else {
                bp1Clipping.add(Clipping.leftClip(pair.getLeft()));
            }
            if (variant.orientationBP2() > 0) {
                bp2Clipping.add(Clipping.rightClip(pair.getRight()));
            } else {
                bp2Clipping.add(Clipping.leftClip(pair.getRight()));
            }
        }

        // NERA: Include more clipping information

        for (final Pair<SAMRecord, SAMRecord> pair : clippedProper) {
            if (stream(pair).allMatch(r -> Location.fromSAMRecord(r).sameChromosomeAs(variant.locationBP1()))) {
                if (variant.orientationBP1() > 0) {
                    bp1Clipping.add(Clipping.rightClip(pair.getRight()));
                } else {
                    bp1Clipping.add(Clipping.leftClip(pair.getLeft()));
                }
            }

            if (stream(pair).allMatch(r -> Location.fromSAMRecord(r).sameChromosomeAs(variant.locationBP2()))) {
                if (variant.orientationBP2() > 0) {
                    bp2Clipping.add(Clipping.rightClip(pair.getRight()));
                } else {
                    bp2Clipping.add(Clipping.leftClip(pair.getLeft()));
                }
            }
        }

        // NERA: Include secondary clipping information
        for (final Pair<SAMRecord, SAMRecord> pair : secondaryPairs) {
            if (stream(pair).allMatch(r -> Location.fromSAMRecord(r).sameChromosomeAs(variant.locationBP1()))) {
                if (variant.orientationBP1() > 0) {
                    bp1Clipping.add(Clipping.rightClip(pair.getRight()));
                } else {
                    bp1Clipping.add(Clipping.leftClip(pair.getLeft()));
                }
            }
            if (stream(pair).allMatch(r -> Location.fromSAMRecord(r).sameChromosomeAs(variant.locationBP2()))) {
                if (variant.orientationBP2() > 0) {
                    bp2Clipping.add(Clipping.rightClip(pair.getRight()));
                } else {
                    bp2Clipping.add(Clipping.leftClip(pair.getLeft()));
                }
            }
        }

        // NERA: Determine candidates based on clipping info
        final List<Location> bp1Candidates = bp1Clipping.getSequences()
                .stream()
                //.filter(c -> c.LongestClipSequence.length() >= 5)
                .map(c -> c.alignment)
                .filter(c -> withinRange(c, variant.locationBP1(), variant.uncertaintyBP1()))
                .collect(Collectors.toList());

        if (bp1Candidates.isEmpty()) {
            interesting.stream()
                    .map(Pair::getLeft)
                    .map(r -> Location.fromSAMRecord(r, variant.orientationBP1() < 0).add(variant.orientationBP1() > 0 ? 1 : -1))
                    .filter(l -> withinRange(l, variant.locationBP1(), variant.uncertaintyBP1()))
                    .max((a, b) -> variant.orientationBP1() > 0 ? a.compareTo(b) : b.compareTo(a))
                    .ifPresent(bp1Candidates::add);
        }

        final List<Location> bp2Candidates = bp2Clipping.getSequences()
                .stream()
                //.filter(c -> c.LongestClipSequence.length() >= 5)
                .map(c -> c.alignment)
                .filter(c -> withinRange(c, variant.locationBP2(), variant.uncertaintyBP2()))
                .collect(Collectors.toList());

        if (bp2Candidates.isEmpty()) {
            interesting.stream()
                    .map(Pair::getRight)
                    .map(r -> Location.fromSAMRecord(r, variant.orientationBP2() < 0).add(variant.orientationBP2() > 0 ? 1 : -1))
                    .filter(l -> withinRange(l, variant.locationBP2(), variant.uncertaintyBP2()))
                    .max((a, b) -> variant.orientationBP2() > 0 ? a.compareTo(b) : b.compareTo(a))
                    .ifPresent(bp2Candidates::add);
        }

        // NERA: NOTE - we include homology on both sides here and take it out later
        final Location breakpoint1 = bp1Candidates.isEmpty() ? null : bp1Candidates.get(0).add(-variant.orientationBP1());
        final Location breakpoint2 = bp2Candidates.isEmpty() ? null : bp2Candidates.get(0).add(-variant.orientationBP2());

        return BreakpointResult.from(Pair.of(breakpoint1, breakpoint2));
    }

    @NotNull
    private static SampleStats collectEvidence(final EnrichedVariantContext variant, final SamReader reader,
            final Pair<Location, Location> breakpoints) {
        final Location bp1 = breakpoints.getLeft();
        final Location bp2 = breakpoints.getRight();

        final SampleStats result = new SampleStats();
        final Pair<Integer, Integer> variantOrientation = Pair.of(variant.orientationBP1(), variant.orientationBP2());

        final List<SAMRecord> currentReads = Lists.newArrayList();
        final SAMRecordIterator iterator = reader.iterator();

        final boolean srOnly = variant.isShortVariant() || variant.isInsert();

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

                final boolean correctOrientation = orientation(pair).equals(variantOrientation);
                final boolean correctChromosome = Location.fromSAMRecord(pair.getLeft()).sameChromosomeAs(variant.locationBP1()) && Location
                        .fromSAMRecord(pair.getRight())
                        .sameChromosomeAs(variant.locationBP2());

                final int MAX_INTRA_PAIR_LENGTH = 400;
                final boolean intraPairLength = (variant.orientationBP1() > 0
                        ? bp1.position() - pair.getLeft().getAlignmentEnd()
                        : pair.getLeft().getAlignmentStart() - bp1.position()) + (variant.orientationBP2() > 0
                        ? breakpoints.getRight().position() - pair.getRight().getAlignmentEnd()
                        : pair.getRight().getAlignmentStart() - bp2.position()) < MAX_INTRA_PAIR_LENGTH;

                boolean isPairEvidence = correctOrientation && correctChromosome && intraPairLength;
                if (isPairEvidence) {
                    final int leftOuter = Location.fromSAMRecord(pair.getLeft(), variant.orientationBP1() > 0).compareTo(bp1);
                    final int rightOuter = Location.fromSAMRecord(pair.getRight(), variant.orientationBP2() > 0).compareTo(bp2);

                    if (variant.orientationBP1() > 0) {
                        isPairEvidence &= leftOuter < 0;
                    } else {
                        isPairEvidence &= leftOuter > 0;
                    }
                    if (variant.orientationBP2() > 0) {
                        isPairEvidence &= rightOuter < 0;
                    } else {
                        isPairEvidence &= rightOuter > 0;
                    }
                }

                if (isPairEvidence) {
                    bp1SRSupport |= exactlyClipsBreakpoint(pair.getLeft(), bp1, variant.orientationBP1());
                    bp2SRSupport |= exactlyClipsBreakpoint(pair.getRight(), bp2, variant.orientationBP2());
                    if (!srOnly) {
                        prSupport = true;
                        result.prEvidence.add(pair);
                    }
                }

                if (proper || secondary) {
                    final boolean clipsBP1 = exactlyClipsBreakpoint(variant.orientationBP1() > 0 ? pair.getRight() : pair.getLeft(),
                            bp1,
                            variant.orientationBP1());
                    final boolean clipsBP2 = exactlyClipsBreakpoint(variant.orientationBP2() > 0 ? pair.getRight() : pair.getLeft(),
                            bp2,
                            variant.orientationBP2());

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
                        result.srEvidence.add(pair);
                    }
                }
            }

            // NERA: Increment read counts
            final boolean srSupport = bp1SRSupport || bp2SRSupport;

            if (srSupport && prSupport) {
                result.bp1Stats.prSrSupport++;
            } else if (bp1SRSupport) {
                result.bp1Stats.srOnlySupport++;
            } else if (prSupport) {
                result.bp1Stats.prOnlySupport++;
            }
            if (bp1PRNormal && bp1SRNormal) {
                result.bp1Stats.prSrNormal++;
            } else if (bp1PRNormal && !srOnly) {
                result.bp1Stats.prOnlyNormal++;
            }

            if (srSupport && prSupport) {
                result.bp2Stats.prSrSupport++;
            } else if (bp2SRSupport) {
                result.bp2Stats.srOnlySupport++;
            } else if (prSupport) {
                result.bp2Stats.prOnlySupport++;
            }
            if (bp2PRNormal && bp2SRNormal) {
                result.bp2Stats.prSrNormal++;
            } else if (bp2PRNormal && !srOnly) {
                result.bp2Stats.prOnlyNormal++;
            }
        }

        iterator.close();
        return result;
    }

    @NotNull
    private static File queryNameSortedBAM(final SamReader reader, final QueryInterval[] intervals, final String name) throws IOException {
        final SAMFileHeader header = reader.getFileHeader().clone();
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);

        final File file = File.createTempFile(name, ".bam");
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, file);

        int intervalLength = 0;
        for (QueryInterval interval :intervals) {
            intervalLength += (1 + interval.end - interval.start);
        }
        double maxReads = intervalLength * MAX_READS_PER_INTERVAL_LENGTH_FOR_DOWNSAMPLING_CUTOFF;

        final List<SAMRecord> records = reader.queryOverlapping(intervals).toList();

        long downsampleIndex = 0;
        if (records.size() > maxReads) {
            downsampleIndex = Math.round((double) records.size() / maxReads);
            if (downsampleIndex > 1) {
                LOGGER.warn(String.format("Downsampling BAM with name %s with index %s", name, downsampleIndex));
            }
        }

        long recordCount = 0;
        long readCount = 0;
        for (SAMRecord record : records) {
            if (downsampleIndex == 0 || (recordCount % downsampleIndex == 0)) {
                readCount++;
                writer.addAlignment(record);
            }
            recordCount++;
        }

        writer.close();
        LOGGER.info(String.format("    Finished writing file. Loaded %s reads from %s records", readCount, recordCount));

        return file;
    }

    private static int orientation(@NotNull final SAMRecord record) {
        return record.getReadNegativeStrandFlag() ? -1 : 1;
    }

    @NotNull
    private static Pair<Integer, Integer> orientation(@NotNull final Pair<SAMRecord, SAMRecord> pair) {
        return Pair.of(orientation(pair.getLeft()), orientation(pair.getRight()));
    }

    @NotNull
    private static <L> Stream<L> stream(@NotNull final Pair<L, L> pair) {
        return Stream.of(pair.getLeft(), pair.getRight());
    }

    private static boolean isMate(@NotNull final SAMRecord read, @NotNull final SAMRecord mate) {
        return read.getReadName().equals(mate.getReadName()) && read.getMateReferenceIndex().equals(mate.getReferenceIndex())
                && Math.abs(read.getMateAlignmentStart() - mate.getAlignmentStart()) <= 1;
    }

    private static boolean span(@NotNull final Pair<SAMRecord, SAMRecord> pair, @NotNull final Location breakpoint) {
        return Location.fromSAMRecord(pair.getLeft(), true).compareTo(breakpoint) <= 0
                && Location.fromSAMRecord(pair.getRight(), false).compareTo(breakpoint) >= 0;
    }

    private static boolean overlap(@NotNull final SAMRecord read, @NotNull final Location breakpoint) {
        return read.getReferenceIndex() == breakpoint.referenceIndex() && read.getAlignmentStart() <= breakpoint.position()
                && breakpoint.position() <= read.getAlignmentEnd();
    }

    private static boolean overlap(@NotNull final Pair<SAMRecord, SAMRecord> pair, @NotNull final Location breakpoint) {
        return stream(pair).anyMatch(r -> overlap(r, breakpoint));
    }

    private static boolean clippedOnCorrectSide(@NotNull final SAMRecord record, final int orientation) {
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

                // NERA: check both directions due to secondary alignments
                if (isMate(r0, r1) || isMate(r1, r0)) {
                    pairs.add(Pair.of(r0, r1));
                }
            }
        }
        return pairs;
    }

    enum BreakpointError {
        NONE,
        ALGO_ERROR
    }

    private static class PairedReads extends ArrayList<Pair<SAMRecord, SAMRecord>> {

    }

    private static class BreakpointResult {

        @NotNull
        private final Pair<Location, Location> breakpoints;
        @NotNull
        private final BreakpointError error;

        private BreakpointResult(@NotNull final Pair<Location, Location> breakpoints) {
            this.breakpoints = breakpoints;
            if (stream(breakpoints).anyMatch(Objects::isNull)) {
                this.error = BreakpointError.ALGO_ERROR;
            } else {
                this.error = BreakpointError.NONE;
            }
        }

        @NotNull
        private static BreakpointResult from(@NotNull final Pair<Location, Location> breakpoints) {
            return new BreakpointResult(breakpoints);
        }
    }
}
