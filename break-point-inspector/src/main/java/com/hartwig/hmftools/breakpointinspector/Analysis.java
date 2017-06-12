package com.hartwig.hmftools.breakpointinspector;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.*;

import static com.hartwig.hmftools.breakpointinspector.ReadHelpers.*;
import static com.hartwig.hmftools.breakpointinspector.Util.*;
import static com.hartwig.hmftools.breakpointinspector.Stats.*;

import org.jetbrains.annotations.Nullable;

class Analysis {

    private static boolean assessDEL(final HMFVariantContext ctx, final List<ReadInfo> pair) {

        final Location adjustedBP1 =
                ctx.Uncertainty1 == null ? ctx.Breakpoint1 : ctx.Breakpoint1.add(ctx.Uncertainty1.End);
        final Location adjustedBP2 =
                ctx.Uncertainty2 == null ? ctx.Breakpoint2 : ctx.Breakpoint2.add(ctx.Uncertainty2.Start);

        boolean pairEvidence;
        // side
        pairEvidence = pair.get(0).Read.getAlignmentStart() <= adjustedBP1.Position;
        pairEvidence &= pair.get(1).Read.getAlignmentEnd() >= adjustedBP2.Position;
        // orientation
        pairEvidence &= !pair.get(0).Read.getReadNegativeStrandFlag() && pair.get(1).Read.getReadNegativeStrandFlag();

        return pairEvidence;
    }

    private static boolean assessDUP(final HMFVariantContext ctx, final List<ReadInfo> pair) {

        final Location adjustedBP1 =
                ctx.Uncertainty1 == null ? ctx.Breakpoint1 : ctx.Breakpoint1.add(ctx.Uncertainty1.Start);
        final Location adjustedBP2 =
                ctx.Uncertainty2 == null ? ctx.Breakpoint2 : ctx.Breakpoint2.add(ctx.Uncertainty2.End);

        boolean pairEvidence;
        // side
        pairEvidence = pair.get(0).Read.getAlignmentEnd() >= adjustedBP1.Position;
        pairEvidence &= pair.get(1).Read.getAlignmentStart() <= adjustedBP2.Position;
        // orientation
        pairEvidence &= pair.get(0).Read.getReadNegativeStrandFlag() && !pair.get(1).Read.getReadNegativeStrandFlag();

        return pairEvidence;
    }

    private static boolean assessINV(final HMFVariantContext ctx, final List<ReadInfo> pair, boolean inv3) {

        boolean pairEvidence;
        if (inv3) {
            final Location adjustedBP1 =
                    ctx.Uncertainty1 == null ? ctx.Breakpoint1 : ctx.Breakpoint1.add(ctx.Uncertainty1.End);
            final Location adjustedBP2 =
                    ctx.Uncertainty2 == null ? ctx.Breakpoint2 : ctx.Breakpoint2.add(ctx.Uncertainty2.End);

            pairEvidence = pair.get(0).Read.getAlignmentStart() <= adjustedBP1.Position
                    && pair.get(1).Read.getAlignmentStart() <= adjustedBP2.Position;
            pairEvidence &= !pair.get(0).Read.getReadNegativeStrandFlag(); // forward strand x ---->
            pairEvidence &= !pair.get(1).Read.getReadNegativeStrandFlag(); // forward strand x ---->
        } else {
            final Location adjustedBP1 =
                    ctx.Uncertainty1 == null ? ctx.Breakpoint1 : ctx.Breakpoint1.add(ctx.Uncertainty1.Start);
            final Location adjustedBP2 =
                    ctx.Uncertainty2 == null ? ctx.Breakpoint2 : ctx.Breakpoint2.add(ctx.Uncertainty2.Start);

            pairEvidence = pair.get(0).Read.getAlignmentEnd() >= adjustedBP1.Position
                    && pair.get(1).Read.getAlignmentEnd() >= adjustedBP2.Position;
            pairEvidence &= pair.get(0).Read.getReadNegativeStrandFlag(); // reverse strand <---- x
            pairEvidence &= pair.get(1).Read.getReadNegativeStrandFlag(); // reverse strand <---- x
        }

        return pairEvidence;
    }

    private static Stats.Sample calculateEvidenceStats(final ClassifiedReadResults queryResult,
            final HMFVariantContext ctx) {

        final Stats.Sample result = new Stats.Sample();

        final Stats.ClipStats pairedClipsAtBP1 = new Stats.ClipStats();
        final Stats.ClipStats pairedClipsAtBP2 = new Stats.ClipStats();

        final List<List<ReadInfo>> spanningPairs = new ArrayList<>();
        final List<List<ReadInfo>> localPairs = new ArrayList<>();

        for (final NamedReadCollection collection : queryResult.ReadMap.values()) {
            for (final ReadInfo p0 : collection.Reads) {

                // update clipping stats
                result.Clipping_Stats.addToClippingStats(p0.Read);

                // find the mate
                final ReadInfo p1 = collection.Reads.stream().filter(i -> isMate(p0.Read, i.Read)).findFirst().orElse(
                        null);

                // single read
                if (p1 == null) {
                    if (p0.Category == ReadCategory.MATE_UNMAPPED) {
                        result.Get(p0.Breakpoint).Unmapped_Mate++;
                    } else if (p0.Category != ReadCategory.NORMAL
                            && p0.Location != Overlap.FILTERED) { // TODO: check condition
                        result.Get(p0.Breakpoint).Diff_Variant++;
                    }
                    continue;
                }

                // special case for BND pairings, we want them in same order as manta breakpoint
                if (p0.Read.getInferredInsertSize() == 0 && !p0.Read.getReferenceIndex().equals(
                        ctx.Breakpoint1.ReferenceIndex)) {
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
                    boolean pairEvidence = false;
                    switch (ctx.Type) {
                        case DEL:
                            pairEvidence = assessDEL(ctx, pair);
                            break;
                        case DUP:
                            pairEvidence = assessDUP(ctx, pair);
                            break;
                        case INV3:
                            pairEvidence = assessINV(ctx, pair, true);
                            break;
                        case INV5:
                            pairEvidence = assessINV(ctx, pair, false);
                            break;
                    }

                    if (pairEvidence) {
                        spanningPairs.add(pair);
                        pairedClipsAtBP1.addToClippingStats(p0.Read);
                        pairedClipsAtBP2.addToClippingStats(p1.Read);
                    } else {
                        result.BP1_Stats.Diff_Variant++;
                        result.BP2_Stats.Diff_Variant++;
                    }
                } else {
                    // p0 == p1
                    localPairs.add(pair);
                }
            }
        }

        // 1: look at clipped pairs first
        if (result.BP1 == null) {
            result.BP1 = pairedClipsAtBP1.LocationMap.entrySet().stream().max(
                    Comparator.comparingInt(a -> a.getValue().Reads.size())).map(Map.Entry::getKey).orElse(null);
        }
        if (result.BP2 == null) {
            result.BP2 = pairedClipsAtBP2.LocationMap.entrySet().stream().max(
                    Comparator.comparingInt(a -> a.getValue().Reads.size())).map(Map.Entry::getKey).orElse(null);
        }

        // 2: then we should look at unpaired clips to determine breakpoint
        if (result.BP1 == null) {
            result.BP1 = result.Clipping_Stats.LocationMap.entrySet().stream().filter(
                    e -> e.getKey().closeTo(ctx.Breakpoint1, ctx.Uncertainty1)).max(
                    Comparator.comparingInt(a -> a.getValue().Reads.size())).map(Map.Entry::getKey).orElse(null);
        }
        if (result.BP2 == null) {
            result.BP2 = result.Clipping_Stats.LocationMap.entrySet().stream().filter(
                    e -> e.getKey().closeTo(ctx.Breakpoint2, ctx.Uncertainty1)).max(
                    Comparator.comparingInt(a -> a.getValue().Reads.size())).map(Map.Entry::getKey).orElse(null);
        }
        // TODO: 3: then we should just use the Manta breakpoint

        for (final List<ReadInfo> pair : spanningPairs) {
            for (final ReadInfo r : pair) {
                final boolean sr = getClips(r.Read).stream().anyMatch(
                        c -> c.Alignment.closeTo(r.Breakpoint == Region.BP1 ? result.BP1 : result.BP2));
                // TODO: or exact equals?
                // TODO: also check side of clip
                if (sr) {
                    result.Get(r.Breakpoint).PR_SR_Support++;
                } else {
                    result.Get(r.Breakpoint).PR_Only_Support++;
                }
            }
        }

        for (final List<ReadInfo> pair : localPairs) {
            final ReadInfo p0 = pair.get(0);
            final ReadInfo p1 = pair.get(1);
            final Stats.BreakPoint stats = result.Get(p0.Breakpoint);

            final boolean sr = pair.stream().anyMatch(r -> getClips(r.Read).stream().anyMatch(
                    c -> c.Alignment.closeTo(r.Breakpoint == Region.BP1 ? result.BP1 : result.BP2)));
            // TODO: also check side of clip
            if (sr) {
                stats.SR_Only_Support++;
            } else if (p0.Location == Overlap.STRADDLE && p1.Location == Overlap.STRADDLE) {
                stats.PR_Only_Normal++;
            } else if (p0.Location == Overlap.INTERSECT || p1.Location == Overlap.INTERSECT) {
                stats.PR_SR_Normal++;
            }
        }

        return result;
    }

    private static Stats.Sample calculateStats(final ClassifiedReadResults queryResult, final HMFVariantContext ctx) {
        return calculateEvidenceStats(queryResult, ctx);
    }

    // TODO: this should really be per sample
    private static Set<Integer> sWrittenReads = new HashSet<>();

    private static ClassifiedReadResults performQueryAndClassify(final SamReader reader,
            @Nullable SAMFileWriter evidenceWriter, final QueryInterval[] intervals, final HMFVariantContext ctx) {

        final ClassifiedReadResults result = new ClassifiedReadResults();

        // execute query and parse the results
        final SAMRecordIterator results = reader.query(intervals, false);
        while (results.hasNext()) {

            final SAMRecord read = results.next();

            // write to evidence bam
            if (evidenceWriter != null) {
                final int hash = read.getSAMString().hashCode();
                if (!sWrittenReads.contains(hash)) {
                    sWrittenReads.add(hash);
                    evidenceWriter.addAlignment(read);
                }
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
            if (!getClips(read).isEmpty()) {
                info.Location = Overlap.CLIP;
            } else {
                info.Location = Overlap.PROXIMITY;
            }
        }
        results.close();

        return result;
    }

    private static boolean compareString(final String a, final String b) {
        if (a.length() < b.length())
            return b.startsWith(a) || b.endsWith(a);
        else
            return a.startsWith(b) || a.endsWith(b);
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

        final List<String> filters = new ArrayList<>(ctx.Filter);
        if (refStats.BP1_Stats.PR_Only_Support > 0 || refStats.BP1_Stats.PR_SR_Support > 0) {
            filters.add("HMF_PRNormalSupport");
        }

        // TODO: check the sequences match too (not just location)
        boolean concordance;
        concordance = refStats.Clipping_Stats.LocationMap.entrySet().stream().anyMatch(
                l -> l.getKey().closeTo(tumorStats.BP1));
        concordance |= refStats.Clipping_Stats.LocationMap.entrySet().stream().anyMatch(
                l -> l.getKey().closeTo(tumorStats.BP2));
        if (concordance)
            filters.add("HMF_ClippingConcordance");

        // output

        final ArrayList<String> data = new ArrayList<>(extraData);
        data.addAll(refStats.GetData());
        data.addAll(tumorStats.GetData());
        data.add(tumorStats.BP1 != null ? tumorStats.BP1.toString() : "err");
        data.add(tumorStats.BP2 != null ? tumorStats.BP2.toString() : "err");
        data.add(filters.isEmpty() ? "PASS" : String.join(";", filters));
        data.add(tumorStats.Clipping_Stats.toString());
        System.out.println(String.join("\t", data));
    }
}
