package com.hartwig.hmftools.breakpointinspector;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import htsjdk.samtools.*;

import static com.hartwig.hmftools.breakpointinspector.ReadHelpers.*;
import static com.hartwig.hmftools.breakpointinspector.Util.*;
import static com.hartwig.hmftools.breakpointinspector.Stats.*;

import org.jetbrains.annotations.Nullable;

class Analysis {

    private static class StructuralVariantContext {
        Location Breakpoint1;
        Location Breakpoint2;
        VariantType Type;

        StructuralVariantContext(final Location bp1, final Location bp2, final VariantType type) {
            Breakpoint1 = bp1;
            Breakpoint2 = bp2;
            Type = type;
        }
    }

    private static Stats.ClipStats calculateClippingStats(final ClassifiedReadResults queryResult) {
        final Stats.ClipStats result = new Stats.ClipStats();
        for (final NamedReadCollection collection : queryResult.ReadMap.values()) {
            for (final ReadInfo info : collection.Reads) {
                final SAMRecord read = info.Read;
                for (final ClipInfo clip : getClips(read)) {

                    final Stats.Clip stats = result.LocationMap.computeIfAbsent(clip.Alignment, k -> new Stats.Clip());

                    if (clip.HardClipped) {
                        stats.HardClippedReads.add(read);
                    } else {
                        if (clip.Sequence.length() > stats.LongestClipSequence.length() && clip.Sequence.contains(
                                stats.LongestClipSequence)) {
                            // the existing sequence supports the new sequence
                            stats.LongestClipSequence = clip.Sequence;
                        } else if (!stats.LongestClipSequence.contains(clip.Sequence)) {
                            // this read does not support the existing sequence
                            continue;
                        }

                        stats.Reads.add(read);
                    }
                }
            }
        }
        return result;
    }

    private static boolean assessDEL(final StructuralVariantContext ctx, final List<ReadInfo> pair,
            final Stats.Sample result) {

        boolean pairEvidence;
        pairEvidence = pair.get(0).Read.getAlignmentStart() <= ctx.Breakpoint1.Position;
        pairEvidence &= pair.get(1).Read.getAlignmentStart() >= ctx.Breakpoint2.Position;
        pairEvidence &= SamPairUtil.getPairOrientation(pair.get(0).Read) == SamPairUtil.PairOrientation.FR;

        for (final ReadInfo r : pair) {
            final boolean clipped = r.Location == Overlap.CLIP_EXACT;
            if (clipped && pairEvidence) {
                result.Get(r.Breakpoint).PR_SR_Support++;
            } else if (pairEvidence) {
                result.Get(r.Breakpoint).PR_Only_Support++;
            } else {
                result.Get(r.Breakpoint).Diff_Variant++;
            }
        }

        return false;
    }

    private static boolean assessDUP(final StructuralVariantContext ctx, final List<ReadInfo> pair,
            final Stats.Sample result) {

        boolean pairEvidence;
        pairEvidence = pair.get(0).Read.getAlignmentStart() >= ctx.Breakpoint1.Position;
        pairEvidence &= pair.get(1).Read.getAlignmentStart() <= ctx.Breakpoint2.Position;
        pairEvidence &= SamPairUtil.getPairOrientation(pair.get(0).Read) == SamPairUtil.PairOrientation.RF;

        for (final ReadInfo r : pair) {
            final boolean clipped = r.Location == Overlap.CLIP_EXACT;
            if (clipped && pairEvidence) {
                result.Get(r.Breakpoint).PR_SR_Support++;
            } else if (pairEvidence) {
                result.Get(r.Breakpoint).PR_Only_Support++;
            } else {
                result.Get(r.Breakpoint).Diff_Variant++;
            }
        }

        return false;
    }

    private static boolean assessINV(final StructuralVariantContext ctx, final List<ReadInfo> pair,
            final Stats.Sample result, boolean inv3) {

        boolean pairEvidence;
        if (inv3) {
            pairEvidence = pair.get(0).Read.getAlignmentStart() <= ctx.Breakpoint1.Position
                    && pair.get(1).Read.getAlignmentStart() <= ctx.Breakpoint2.Position;
        } else {
            pairEvidence = pair.get(0).Read.getAlignmentStart() >= ctx.Breakpoint1.Position
                    && pair.get(1).Read.getAlignmentStart() >= ctx.Breakpoint2.Position;
        }
        // TODO: determine tandem direction (L-L or R-R)
        pairEvidence &= SamPairUtil.getPairOrientation(pair.get(0).Read) == SamPairUtil.PairOrientation.TANDEM;

        for (final ReadInfo r : pair) {
            final boolean clipped = r.Location == Overlap.CLIP_EXACT;
            if (clipped && pairEvidence) {
                result.Get(r.Breakpoint).PR_SR_Support++;
            } else if (pairEvidence) {
                result.Get(r.Breakpoint).PR_Only_Support++;
            } else {
                result.Get(r.Breakpoint).Diff_Variant++;
            }
        }

        return false;
    }

    private static Stats.Sample calculateEvidenceStats(final ClassifiedReadResults queryResult,
            final StructuralVariantContext ctx) {
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

                // at this stage, we have a legitimate pair
                final List<ReadInfo> pair = Arrays.asList(p0, p1);

                // supports the break point
                if (p0.Breakpoint != p1.Breakpoint) {
                    switch (ctx.Type) {
                        case DEL:
                            assessDEL(ctx, pair, result);
                            break;
                        case DUP:
                            assessDUP(ctx, pair, result);
                            break;
                        case INV3:
                            assessINV(ctx, pair, result, true);
                            break;
                        case INV5:
                            assessINV(ctx, pair, result, false);
                            break;
                    }
                } else {
                    final Stats.BreakPoint stats = result.Get(p0.Breakpoint);
                    final boolean splitEvidence = pair.stream().anyMatch(p -> p.Location == Overlap.CLIP_EXACT);
                    if (splitEvidence) {
                        stats.SR_Only_Support++;
                    } else if (p0.Location == Overlap.STRADDLE && p1.Location == Overlap.STRADDLE) {
                        stats.PR_Only_Normal++;
                    } else if (p0.Location == Overlap.INTERSECT || p1.Location == Overlap.INTERSECT) {
                        stats.PR_SR_Normal++;
                    }
                }
            }
        }
        return result;
    }

    private static Stats.Sample calculateStats(final ClassifiedReadResults queryResult,
            final StructuralVariantContext ctx) {
        final Stats.Sample result = calculateEvidenceStats(queryResult, ctx);
        result.Clipping_Stats = calculateClippingStats(queryResult);
        return result;
    }

    private static ClassifiedReadResults performQueryAndClassify(final SamReader reader,
            @Nullable SAMFileWriter evidenceWriter, final QueryInterval[] intervals,
            final StructuralVariantContext ctx) {

        final ClassifiedReadResults result = new ClassifiedReadResults();

        // execute query and parse the results
        final SAMRecordIterator results = reader.query(intervals, false);
        while (results.hasNext()) {

            final SAMRecord read = results.next();
            if (evidenceWriter != null) {
                evidenceWriter.addAlignment(read);
            }

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

            info.Breakpoint = determineRegion(read, ctx.Breakpoint1, ctx.Breakpoint2);
            if (info.Breakpoint == Region.OTHER) {
                throw new RuntimeException("read from unexpected region");
            }

            final Location bp = info.Breakpoint == Region.BP1 ? ctx.Breakpoint1 : ctx.Breakpoint2;
            final boolean normal = read.getProperPairFlag();
            final boolean clipped = isClipped(info.Read);

            if (normal) {
                info.Category = ReadCategory.NORMAL;

                if (!clipped) {
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
            final ClipInfo left = getLeftClip(read);
            final ClipInfo right = getRightClip(read);

            // TODO: double check clipping conditions
            if (left != null && left.Alignment.add(-1).equals(bp)) {
                info.Location = Overlap.CLIP_EXACT;
            } else if (right != null && right.Alignment.equals(bp)) {
                info.Location = Overlap.CLIP_EXACT;
            } else if (left != null || right != null) {
                info.Location = Overlap.CLIP_OTHER;
            } else {
                info.Location = Overlap.PROXIMITY;
            }
        }
        results.close();

        return result;
    }

    static void processStructuralVariant(final List<String> extraData, final SamReader refReader,
            @Nullable final SAMFileWriter refWriter, final SamReader tumorReader,
            @Nullable final SAMFileWriter tumorWriter, final Location location1, final Location location2,
            final int range, VariantType svType) {
        // work out the query intervals
        QueryInterval[] queryIntervals = {
                new QueryInterval(location1.ReferenceIndex, Math.max(0, location1.Position - range),
                        location1.Position + range),
                new QueryInterval(location2.ReferenceIndex, Math.max(0, location2.Position - range),
                        location2.Position + range) };
        queryIntervals = QueryInterval.optimizeIntervals(queryIntervals);

        // begin processing
        final StructuralVariantContext context = new StructuralVariantContext(location1, location2, svType);

        final ClassifiedReadResults refResult = performQueryAndClassify(refReader, refWriter, queryIntervals, context);
        final Sample refStats = calculateStats(refResult, context);

        final ClassifiedReadResults tumorResult = performQueryAndClassify(tumorReader, tumorWriter, queryIntervals,
                context);
        final Sample tumorStats = calculateStats(tumorResult, context);

        final ArrayList<String> data = new ArrayList<>(extraData);
        data.addAll(refStats.GetData());
        data.addAll(tumorStats.GetData());
        data.add(tumorStats.Clipping_Stats.toString());
        System.out.println(String.join("\t", data));
    }
}
