package com.hartwig.hmftools.breakpointinspector;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Objects;

import htsjdk.samtools.*;

import static com.hartwig.hmftools.breakpointinspector.ReadHelpers.*;
import static com.hartwig.hmftools.breakpointinspector.Util.*;
import static com.hartwig.hmftools.breakpointinspector.Stats.*;

import org.jetbrains.annotations.Nullable;

class Analysis {

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

    private static boolean checkSR(final HMFVariantContext ctx, final ReadInfo r) {
        final Location breakpoint = r.Breakpoint == Region.BP1 ? ctx.Breakpoint1 : ctx.Breakpoint2;
        final Range uncertainty = r.Breakpoint == Region.BP1 ? ctx.Uncertainty1 : ctx.Uncertainty2;

        if (r.Location == Overlap.CLIP_EXACT) {
            return true;
        } else if (r.Location == Overlap.CLIP_OTHER && uncertainty != null) {
            final ClipInfo left = getLeftClip(r.Read);
            if (left != null && left.Alignment.Position >= (breakpoint.Position + uncertainty.Min)
                    && left.Alignment.Position <= (breakpoint.Position + uncertainty.Max)) {
                return true;
            }
            final ClipInfo right = getRightClip(r.Read);
            if (right != null && right.Alignment.Position >= (breakpoint.Position + uncertainty.Min)
                    && right.Alignment.Position <= (breakpoint.Position + uncertainty.Max)) {
                return true;
            }
        }

        return false;
    }

    private static void assessDEL(final HMFVariantContext ctx, final List<ReadInfo> pair, final Stats.Sample result) {

        boolean pairEvidence;
        pairEvidence = pair.get(0).Read.getAlignmentStart() <= ctx.Breakpoint1.Position;
        pairEvidence &= pair.get(1).Read.getAlignmentEnd() >= ctx.Breakpoint2.Position;

        if (Objects.equals(pair.get(0).Read.getReferenceIndex(), pair.get(1).Read.getReferenceIndex())) {
            // same chromosome
            pairEvidence &= SamPairUtil.getPairOrientation(pair.get(0).Read) == SamPairUtil.PairOrientation.FR;
        } else {
            // handle BND case
            pairEvidence &=
                    pair.get(0).Read.getReadNegativeStrandFlag() != pair.get(1).Read.getReadNegativeStrandFlag();
        }

        for (final ReadInfo r : pair) {
            if (checkSR(ctx, r) && pairEvidence) {
                result.Get(r.Breakpoint).PR_SR_Support++;
            } else if (pairEvidence) {
                result.Get(r.Breakpoint).PR_Only_Support++;
            } else {
                result.Get(r.Breakpoint).Diff_Variant++;
            }
        }
    }

    private static void assessDUP(final HMFVariantContext ctx, final List<ReadInfo> pair, final Stats.Sample result) {

        boolean pairEvidence;
        pairEvidence = pair.get(0).Read.getAlignmentEnd() >= ctx.Breakpoint1.Position;
        pairEvidence &= pair.get(1).Read.getAlignmentStart() <= ctx.Breakpoint2.Position;
        pairEvidence &= SamPairUtil.getPairOrientation(pair.get(0).Read) == SamPairUtil.PairOrientation.RF;

        for (final ReadInfo r : pair) {
            if (checkSR(ctx, r) && pairEvidence) {
                result.Get(r.Breakpoint).PR_SR_Support++;
            } else if (pairEvidence) {
                result.Get(r.Breakpoint).PR_Only_Support++;
            } else {
                result.Get(r.Breakpoint).Diff_Variant++;
            }
        }
    }

    private static void assessINV(final HMFVariantContext ctx, final List<ReadInfo> pair, final Stats.Sample result,
            boolean inv3) {

        boolean pairEvidence;
        if (inv3) {
            pairEvidence = pair.get(0).Read.getAlignmentStart() <= ctx.Breakpoint1.Position
                    && pair.get(1).Read.getAlignmentStart() <= ctx.Breakpoint2.Position;
            pairEvidence &= !pair.get(0).Read.getReadNegativeStrandFlag(); // forward strand x ---->
            pairEvidence &= !pair.get(1).Read.getReadNegativeStrandFlag(); // forward strand x ---->
        } else {
            pairEvidence = pair.get(0).Read.getAlignmentEnd() >= ctx.Breakpoint1.Position
                    && pair.get(1).Read.getAlignmentEnd() >= ctx.Breakpoint2.Position;
            pairEvidence &= pair.get(0).Read.getReadNegativeStrandFlag(); // reverse strand <---- x
            pairEvidence &= pair.get(1).Read.getReadNegativeStrandFlag(); // reverse strand <---- x
        }

        for (final ReadInfo r : pair) {
            if (checkSR(ctx, r) && pairEvidence) {
                result.Get(r.Breakpoint).PR_SR_Support++;
            } else if (pairEvidence) {
                result.Get(r.Breakpoint).PR_Only_Support++;
            } else {
                result.Get(r.Breakpoint).Diff_Variant++;
            }
        }
    }

    private static Stats.Sample calculateEvidenceStats(final ClassifiedReadResults queryResult,
            final HMFVariantContext ctx) {
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
                    final boolean splitEvidence = pair.stream().anyMatch(p -> checkSR(ctx, p));
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

    private static Stats.Sample calculateStats(final ClassifiedReadResults queryResult, final HMFVariantContext ctx) {
        final Stats.Sample result = calculateEvidenceStats(queryResult, ctx);
        result.Clipping_Stats = calculateClippingStats(queryResult);
        return result;
    }

    private static ClassifiedReadResults performQueryAndClassify(final SamReader reader,
            @Nullable SAMFileWriter evidenceWriter, final QueryInterval[] intervals, final HMFVariantContext ctx) {

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
            if (left != null && left.Alignment.closeTo(bp)) {
                info.Location = Overlap.CLIP_EXACT;
            } else if (right != null && right.Alignment.closeTo(bp)) {
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
            @Nullable final SAMFileWriter tumorWriter, final HMFVariantContext ctx, final int range) {
        // work out the query intervals
        QueryInterval[] queryIntervals = {
                new QueryInterval(ctx.Breakpoint1.ReferenceIndex, Math.max(0, ctx.Breakpoint1.Position - range),
                        ctx.Breakpoint1.Position + range),
                new QueryInterval(ctx.Breakpoint2.ReferenceIndex, Math.max(0, ctx.Breakpoint2.Position - range),
                        ctx.Breakpoint2.Position + range) };
        queryIntervals = QueryInterval.optimizeIntervals(queryIntervals);

        // begin processing

        final ClassifiedReadResults refResult = performQueryAndClassify(refReader, refWriter, queryIntervals, ctx);
        final Sample refStats = calculateStats(refResult, ctx);

        final ClassifiedReadResults tumorResult = performQueryAndClassify(tumorReader, tumorWriter, queryIntervals,
                ctx);
        final Sample tumorStats = calculateStats(tumorResult, ctx);

        // filtering

        String filterReason = "PASS";
        if (refStats.BP1_Stats.PR_Only_Support > 0 || refStats.BP1_Stats.PR_SR_Support > 0) {
            filterReason = "NormalSupport";
        }

        for (final Map.Entry<Location, Clip> entry : refStats.Clipping_Stats.LocationMap.entrySet()) {
            final Clip tumorClip = tumorStats.Clipping_Stats.LocationMap.get(entry.getKey());
            if (tumorClip != null) {
                // TODO: if the clip is a supporting read ??
            }
        }

        // output

        final ArrayList<String> data = new ArrayList<>(extraData);
        data.addAll(refStats.GetData());
        data.addAll(tumorStats.GetData());
        data.add(filterReason);
        data.add(tumorStats.Clipping_Stats.toString());
        System.out.println(String.join("\t", data));
    }
}
