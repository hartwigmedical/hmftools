package com.hartwig.hmftools.breakpointinspector;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import htsjdk.samtools.*;

import static com.hartwig.hmftools.breakpointinspector.Util.*;
import static com.hartwig.hmftools.breakpointinspector.Stats.*;

public class Analysis {

    private static final int MANTA_REQ_PAIR_MIN = 50;
    private static final int MANTA_REQ_SPLIT_MIN = 15;

    private static boolean readIntersectsLocation(final SAMRecord read, final Location location) {
        assert !read.getReadUnmappedFlag();
        return read.getReferenceIndex() == location.ReferenceIndex && read.getAlignmentStart() <= location.Position
                && read.getAlignmentEnd() >= location.Position;
    }

    private static boolean pairStraddlesLocation(final SAMRecord read, final Location location) {
        assert !read.getReadUnmappedFlag() && !read.getMateUnmappedFlag();
        if (read.getInferredInsertSize() > 0)
            return read.getAlignmentStart() <= location.Position
                    && (read.getAlignmentStart() + read.getInferredInsertSize()) >= location.Position;
        else
            return read.getMateAlignmentStart() <= location.Position
                    && (read.getMateAlignmentStart() - read.getInferredInsertSize()) >= location.Position;
    }

    private static ClipInfo getClipInfo(final SAMRecord read) {
        final ClipInfo result = new ClipInfo();
        final Cigar cigar = read.getCigar();
        switch (cigar.getFirstCigarElement().getOperator()) {
            case H:
                result.HardClipped = true;
            case S:
                result.Side = ClipSide.LEFT_CLIP;
                result.Length = cigar.getFirstCigarElement().getLength();
                return result;
        }
        switch (cigar.getLastCigarElement().getOperator()) {
            case H:
                result.HardClipped = true;
            case S:
                result.Side = ClipSide.RIGHT_CLIP;
                result.Length = cigar.getLastCigarElement().getLength();
                return result;
        }
        return result;
    }

    private static Region determineRegion(final SAMRecord read, final Location location1, final Location location2) {
        if (location1.ReferenceIndex != location2.ReferenceIndex) {
            if (read.getReferenceIndex() == location1.ReferenceIndex)
                return Region.BP1;
            else if (read.getReferenceIndex() == location2.ReferenceIndex)
                return Region.BP2;
        } else if (read.getReferenceIndex() == location1.ReferenceIndex) {
            if (Math.abs(read.getAlignmentStart() - location1.Position) < Math.abs(
                    read.getAlignmentStart() - location2.Position))
                return Region.BP1;
            else
                return Region.BP2;
        }
        return Region.OTHER;
    }

    private static boolean isMate(final SAMRecord read, final SAMRecord mate) {
        return read.getReadName().equals(mate.getReadName())
                && read.getMateReferenceIndex() == mate.getReferenceIndex()
                && read.getMateAlignmentStart() == mate.getAlignmentStart();
    }

    private static Stats.ClipStats calculateClippingStats(final ClassifiedReadResults queryResult) {
        final Stats.ClipStats result = new Stats.ClipStats();
        for (final NamedReadCollection collection : queryResult.ReadMap.values()) {
            for (final ReadInfo info : collection.Reads) {
                final SAMRecord read = info.Read;
                if (info.Location == Overlap.CLIP) {
                    // TODO: handle multiple clip sides
                    final Location alignment = Location.fromSAMRecord(read, info.Clipping.Side == ClipSide.LEFT_CLIP);
                    final Stats.Clip clip = result.LocationMap.computeIfAbsent(alignment, k -> new Stats.Clip() {{
                        Side = info.Clipping.Side;
                    }});
                    if (info.Clipping.HardClipped) {
                        clip.HardClippedReads.add(read);
                    } else {
                        final String clippedSequence;
                        if (info.Clipping.Side == ClipSide.LEFT_CLIP) {
                            clippedSequence = read.getReadString().substring(0, info.Clipping.Length);
                        } else {
                            clippedSequence = read.getReadString().substring(
                                    read.getReadLength() - info.Clipping.Length);
                        }

                        if (clippedSequence.length() > clip.LongestClipSequence.length() && clippedSequence.contains(
                                clip.LongestClipSequence)) {
                            // the existing sequence supports the new sequence
                            clip.LongestClipSequence = clippedSequence;
                        } else if (!clip.LongestClipSequence.contains(clippedSequence)) {
                            // this read does not support the existing sequence
                            continue;
                        }

                        clip.Reads.add(read);
                    }
                }
            }
        }
        return result;
    }

    private static Stats.Sample calculateEvidenceStats(final ClassifiedReadResults queryResult) {
        final Stats.Sample result = new Stats.Sample();
        for (final NamedReadCollection collection : queryResult.ReadMap.values()) {
            // consider the pairings
            for (final ReadInfo p0 : collection.Reads) {
                // find the mate
                final ReadInfo p1 = collection.Reads.stream().filter(i -> isMate(p0.Read, i.Read)).findFirst().orElse(
                        null);

                // single read
                if (p1 == null) {
                    if (p0.Category == ReadCategory.MATE_UNMAPPED) {
                        result.Get(p0.Breakpoint).Unmapped_Mate++;
                    } else if (p0.Category != ReadCategory.NORMAL || p0.Location != Overlap.FILTERED) {
                        result.Get(p0.Breakpoint).Diff_Variant++;
                    }
                    continue;
                }

                // don't consider pairs twice from the reverse pairing
                if (p0.Read.getInferredInsertSize() < 0) {
                    continue;
                }

                // possible two paired but unmapped reads?
                if (p0.Breakpoint == Region.OTHER && p1.Breakpoint == Region.OTHER) {
                    continue;
                } else if (p1.Breakpoint == Region.OTHER) { // must be unmapped mate
                    result.Get(p0.Breakpoint).Unmapped_Mate++;
                    continue;
                } else if (p0.Breakpoint == Region.OTHER) { // must be unmapped mate
                    result.Get(p1.Breakpoint).Unmapped_Mate++;
                    continue;
                }

                final List<ReadInfo> pair = Arrays.asList(p0, p1);
                final boolean pairPossible = pair.stream().allMatch(p -> p.Clipping.Length <= MANTA_REQ_PAIR_MIN);

                if (p0.Breakpoint != p1.Breakpoint) {
                    // supports the break point
                    boolean interesting = false;
                    for (final ReadInfo r : pair) {
                        final boolean splitEvidence =
                                r.Location == Overlap.CLIP && r.Clipping.Length >= MANTA_REQ_SPLIT_MIN;
                        if (pairPossible && splitEvidence) {
                            result.Get(r.Breakpoint).PR_SR_Support++;
                        } else if (pairPossible) {
                            result.Get(r.Breakpoint).PR_Only_Support++;
                        } else if (splitEvidence) {
                            result.Get(r.Breakpoint).SR_Only_Support++;
                        }

                        interesting |= pairPossible || splitEvidence;
                    }

                    if (interesting) {
                        switch (SamPairUtil.getPairOrientation(p0.Read)) {
                            case FR:
                                result.Orientation_Stats.InnieCount++;
                                break;
                            case RF:
                                result.Orientation_Stats.OutieCount++;
                                break;
                            case TANDEM:
                                result.Orientation_Stats.TandemCount++;
                                break;
                        }
                    }
                } else {
                    final Stats.BreakPoint stats = p0.Breakpoint == Region.BP1 ? result.BP1_Stats : result.BP2_Stats;
                    final boolean splitEvidence = pair.stream().anyMatch(
                            p -> p.Location == Overlap.CLIP && p.Clipping.Length >= MANTA_REQ_SPLIT_MIN);
                    if (splitEvidence) {
                        stats.SR_Only_Support++;
                    } else if (p0.Location == Overlap.STRADDLE && p1.Location == Overlap.STRADDLE) {
                        stats.PR_Only_Normal++;
                    } else if (p0.Location == Overlap.INTERSECT || p1.Location == Overlap.INTERSECT) {
                        stats.PR_SR_Normal++;
                        // TODO: does the intersection have to occur within a certain bound of read?
                    }
                }
            }
        }
        return result;
    }

    private static Stats.Sample calculateStats(final ClassifiedReadResults queryResult) {
        final Stats.Sample result = calculateEvidenceStats(queryResult);
        result.Clipping_Stats = calculateClippingStats(queryResult);
        return result;
    }

    private static ClassifiedReadResults performQueryAndClassify(final SamReader reader,
            final QueryInterval[] intervals, final Location bp1, final Location bp2) {

        final ClassifiedReadResults result = new ClassifiedReadResults();

        // execute query and parse the results
        final SAMRecordIterator results = reader.query(intervals, false);
        while (results.hasNext()) {

            final SAMRecord read = results.next();
            final NamedReadCollection collection = result.ReadMap.computeIfAbsent(read.getReadName(),
                    k -> new NamedReadCollection());
            final ReadInfo info = new ReadInfo();

            info.Read = read;
            collection.Reads.add(info);

            // if unmapped there's nothing to do
            if (read.getReadUnmappedFlag()) {
                info.Breakpoint = Region.OTHER;
                info.Category = ReadCategory.UNMAPPED;
                info.Location = Overlap.FILTERED;
                continue;
            }

            info.Clipping = getClipInfo(read);
            info.Breakpoint = determineRegion(read, bp1, bp2);
            if (info.Breakpoint == Region.OTHER) {
                throw new RuntimeException("read from unexpected region");
            }

            final Location bp = info.Breakpoint == Region.BP1 ? bp1 : bp2;
            final boolean normal = read.getProperPairFlag();

            if (normal) {
                info.Category = ReadCategory.NORMAL;

                if (info.Clipping.Side == ClipSide.NONE) {
                    if (readIntersectsLocation(read, bp)) {
                        info.Location = Overlap.INTERSECT;
                    } else if (pairStraddlesLocation(read, bp)) {
                        info.Location = Overlap.STRADDLE;
                    } else {
                        info.Location = Overlap.FILTERED;
                    }
                    continue;
                }

            } else {
                // determine type
                if (read.getMateUnmappedFlag()) {
                    info.Category = ReadCategory.MATE_UNMAPPED;
                } else if (read.getSupplementaryAlignmentFlag()) {
                    info.Category = ReadCategory.CHIMERIC;
                } else if (read.getNotPrimaryAlignmentFlag()) {
                    info.Category = ReadCategory.SECONDARY;
                } else if (read.getReadPairedFlag())
                    info.Category = ReadCategory.SPAN;
            }

            // determine classification
            // TODO: double check clipping conditions
            if (info.Clipping.Side == ClipSide.LEFT_CLIP && read.getAlignmentStart() - 1 == bp.Position) {
                info.Location = Overlap.CLIP;
            } else if (info.Clipping.Side == ClipSide.RIGHT_CLIP && read.getAlignmentEnd() == bp.Position) {
                info.Location = Overlap.CLIP;
            } else {
                info.Location = Overlap.PROXIMITY;
            }

            // if normal read and we don't clip, then it's not interesting
            if (normal && info.Location == Overlap.PROXIMITY)
                info.Location = Overlap.FILTERED;
        }
        results.close();

        return result;
    }

    static void processStructuralVariant(final List<String> extraData, final SamReader refReader,
            final SamReader tumorReader, final Location location1, final Location location2, final int range) {
        // work out the query intervals
        QueryInterval[] queryIntervals = {
                new QueryInterval(location1.ReferenceIndex, Math.max(0, location1.Position - range),
                        location1.Position + range),
                new QueryInterval(location2.ReferenceIndex, Math.max(0, location2.Position - range),
                        location2.Position + range) };
        queryIntervals = QueryInterval.optimizeIntervals(queryIntervals);

        // begin processing

        final ClassifiedReadResults refResult = performQueryAndClassify(refReader, queryIntervals, location1,
                location2);
        final Sample refStats = calculateStats(refResult);

        final ClassifiedReadResults tumorResult = performQueryAndClassify(tumorReader, queryIntervals, location1,
                location2);
        final Sample tumorStats = calculateStats(tumorResult);

        final ArrayList<String> data = new ArrayList<>(extraData);
        data.addAll(refStats.GetData());
        data.addAll(tumorStats.GetData());
        data.add(tumorStats.Clipping_Stats.toString());
        System.out.println(String.join("\t", data));
    }
}
